#include <Arduino.h>
#include <vector>
#include <array>
#include <ArduinoEigenDense.h>

struct RCACFLAG {
    int RegZ = 1;
    float R0 = 10000;
    float Rz = 1;
    float Ru = 0;
    float FeedForward = 0;
    float lambda = 1;
    int Nc = 3; // 0 based indexing
    int lx = 2; // 1 based indexing, number of states
    int lz = 2; // 1 based indexing, number of errors
    int ly = 2; // 1 based indexing, nummber of output
    int lu = 2; // 1 based indexing, number of inputs
    int Integrator = 0;
    String ControllerType = "PID";
};

struct RCACFILT {
    int Nf = 3;
    float Nu = 1;
    float Du = 0;
};

static RCACFLAG FLAG;
static RCACFILT FILT;

static int ltheta;  // To store the result of CalculateRegSize
static float** u_h; // Control Buffer
static float** z_h; // Performance Buffer
static float** r_h; // Command Buffer
static float** yp_h; // Regressor Buffer
static float* intg; // Integration Buffer (size lz x 1)

// filtering buffer
std::vector<std::vector<std::vector<float>>> PHI_window; // [lu][ltheta][pc+pn]
std::vector<std::vector<float>> u_window;                // [lu][pn]
std::vector<std::vector<float>> z_window;                // [lz][pn]
std::vector<std::vector<std::vector<float>>> PHI_filt_window; // [lz][ltheta][Nf]
std::vector<std::vector<std::vector<float>>> PHI_window_collapse_input;
std::vector<std::vector<float>> u_filt_window;           // [lz][Nf]
std::vector<float> theta_k;                              // [ltheta]

std::vector<std::vector<std::vector<float>>> PHI_b_rr; // [lu][ltheta][2*ltheta]
std::vector<std::vector<float>> U_b_rr;
std::vector<std::vector<float>> u_window_vec_input;

int nf_end = 4; // how far it stores data in the buffer, 0 based indexing
int pc = nf_end;
int pn = pc + FILT.Nf + FLAG.Nc;
Eigen::MatrixXf P_k;
//

int kk = 1;



// these would be defined by the user at runtime, just initalizing them
float u_in[2] = {1,2}; // size lu
float z_in[2] = {2,3}; // size lz
float yp_in[2] = {4,5}; // size lz or ly based on regz
float r_in[2] = {6,7}; // size lz
//


std::vector<std::vector<float>> BuildRegressor() {
    const int lz = FLAG.lz;
    const int ly = FLAG.ly;

    // Define PHI as a dynamic array to handle multiple rows and columns
    std::vector<std::vector<float>> PHI;

    if (FLAG.ControllerType == "PID") {
        // PHI includes current `yp_h` (all rows), integration state, and the difference between consecutive rows of `yp_h`
        for (int row = 0; row < lz; row++) {
            PHI.push_back({yp_h[row][0], intg[row], yp_h[row][0] - yp_h[row][1]});
        }
    } else if (FLAG.ControllerType == "PI") {
        // PHI includes current `yp_h` and integration state
        for (int row = 0; row < lz; row++) {
            PHI.push_back({yp_h[row][0], intg[row]});
        }
    } else if (FLAG.ControllerType == "P") {
        // PHI includes only current `yp_h`
        for (int row = 0; row < lz; row++) {
            PHI.push_back({yp_h[row][0]});
        }
    } else if (FLAG.ControllerType == "FSF") {
        // PHI includes only the current `yp_h` (or state for FSF)
        for (int row = 0; row < ly; row++) {
            PHI.push_back({yp_h[row][0]});
        }
    } else if (FLAG.ControllerType == "FSFI") {
        // PHI includes the current `yp_h` and negative integration state
        for (int row = 0; row < lz; row++) {
            PHI.push_back({yp_h[row][0], -intg[row]});
        }
    }

    // Debugging: Print PHI
    Serial.println("PHI Matrix:");
    for (size_t i = 0; i < PHI.size(); i++) {
        Serial.print("Row ");
        Serial.print(i);
        Serial.print(": ");
        for (size_t j = 0; j < PHI[i].size(); j++) {
            Serial.print(PHI[i][j]);
            Serial.print(" ");
        }
        Serial.println();
    }

    return PHI;
} // Tested and works


void updateWindowBuffer(std::vector<std::vector<float>> PHI){
  int nf = FILT.Nf;
  int lu = FLAG.lu;
  int length_window = nf + nf_end;
  for (int row = 0; row < lu; row++) {
    for (int column = length_window; column > 0; column--) {
        u_window[row][column] = u_window[row][column - 1];
    }
    u_window[row][0] = u_in[row];
  }

  int lz = FLAG.lz;
  for (int row = 0; row < lz; row++) {
    for (int column = length_window; column > 0; column--) {
        z_window[row][column] = z_window[row][column - 1];
    }
    z_window[row][0] = z_in[row];
  }


  for (int row1 = 0; row1 < lz; row1++) {
    for (int row2 = 0; row2 < ltheta; row2++) {
      for (int column = length_window; column > 0; column--) {
        PHI_window[row1][row2][column] = PHI_window[row1][row2][column - 1];
      }
      PHI_window[row1][row2][0] = PHI[row1][row2];
    }
  }
} 

std::vector<std::vector<float>> collapse(const std::vector<std::vector<std::vector<float>>>& A) {
    // Initialize the resulting 2D vector with m*k rows and n columns
    std::vector<std::vector<float>> Acol;
    for (const auto &layer : A) {
        for (const auto &row : layer) {
            Acol.push_back(row);
        }
    }
    return Acol;
}

// Function to convert a 2D matrix to a column vector
std::vector<float> vec(const std::vector<std::vector<float>>& matrix) {
    std::vector<float> columnVector;

    // Iterate through each row and column of the matrix
    for (const auto& row : matrix) {
        for (const float& element : row) {
            columnVector.push_back(element);
        }
    }

    return columnVector;
}

// Function to multiply a matrix by a scalar
std::vector<std::vector<float>> multiplyMatrixByScalar(const std::vector<std::vector<float>>& matrix, float scalar) {
    std::vector<std::vector<float>> result = matrix; // Create a copy of the matrix

    // Iterate through each row and column, multiplying by the scalar
    for (auto& row : result) {
        for (auto& element : row) {
            element *= scalar;
        }
    }

    return result;
}

// Function to subtract two matrices of the same size
std::vector<std::vector<float>> subtractMatrices(const std::vector<std::vector<float>>& matrix1, const std::vector<std::vector<float>>& matrix2) {
    std::vector<std::vector<float>> result = matrix1; // Create a copy of the first matrix

    // Iterate through each row and column, subtracting elements of the second matrix
    for (size_t i = 0; i < matrix1.size(); i++) {
        for (size_t j = 0; j < matrix1[i].size(); j++) {
            result[i][j] -= matrix2[i][j];
        }
    }

    return result;
}

std::vector<std::vector<float>> transposeRowVector(const std::vector<float>& rowVector) {
    std::vector<std::vector<float>> columnVector;
    for (const float& element : rowVector) {
        columnVector.push_back({element}); // Each element becomes a row in the column vector
    }
    return columnVector;
}



void FilterSignals(){
  int m = PHI_window.size();              // Number of rows in each 2D matrix
  int n = PHI_window[0].size();           // Number of columns in each 2D matrix
  int row1_u_window = u_window.size();
  // int k = PHI_window[0][0].size();        // Number of 2D matrices in the tensor
  if(kk == 1){ // resize and intialize to zeros
    PHI_window_collapse_input.resize(m, std::vector<std::vector<float>>(n, std::vector<float>(FILT.Nf, 0.0f)));
    PHI_filt_window.resize(FLAG.lz, std::vector<std::vector<float>>(ltheta, std::vector<float>(FILT.Nf, 0.0f)));
    u_window_vec_input.resize(row1_u_window,std::vector<float>(FILT.Nf, 0.0f));
    u_filt_window.resize(FLAG.lz,std::vector<float>(FILT.Nf, 0.0f));
  }
    for (int row1 = 0; row1 < m; row1++) {
      for (int row2 = 0; row2 < n; row2++) { 
        for (int row3 = 0; row3 < FILT.Nf; row3++)
          PHI_window_collapse_input[row1][row2][row3] = PHI_window[row1][row2][row3+1];
      }
    } 

  std::vector<std::vector<float>> Phi_b_rr = collapse(PHI_window_collapse_input);

  for(int row1 = 0; row1<row1_u_window; row1++){
    for(int row2 = 0; row2<FILT.Nf; row2++){
      u_window_vec_input[row1][row2] = u_window[row1][row2];
    }
  }

  std::vector<std::vector<float>> U_b_rr = transposeRowVector(vec(u_window_vec_input));
  Print3DBuffer("PHI_filt_window",PHI_filt_window);
  std::vector<std::vector<float>> Phi_f_b_rr = collapse(PHI_filt_window);
  Print2DBuffer("Phi_f_b_rr",Phi_f_b_rr);

  Print2DBuffer("u_filt_window",u_filt_window);
  std::vector<std::vector<float>> U_b_f_rr = transposeRowVector(vec(u_filt_window));
  Print2DBuffer("U_b_f_rr",U_b_f_rr);

  std::vector<std::vector<float>> PHI_filt_part1 = multiplyMatrixByScalar(Phi_b_rr, FILT.Nu);
  std::vector<std::vector<float>> PHI_filt_part2 = multiplyMatrixByScalar(Phi_f_b_rr, FILT.Du);
  
  std::vector<std::vector<float>> PHI_filt = subtractMatrices(PHI_filt_part1, PHI_filt_part2);

  std::vector<std::vector<float>> u_filt_part1 = multiplyMatrixByScalar(U_b_rr, FILT.Nu);
  std::vector<std::vector<float>> u_filt_part2 = multiplyMatrixByScalar(U_b_f_rr, FILT.Du);
  
  std::vector<std::vector<float>> u_filt = subtractMatrices(u_filt_part1, u_filt_part2);

  for (int row1 = 0; row1 < FLAG.lz; row1++) {
    for (int row2 = 0; row2 < ltheta; row2++) {
      for (int column = FILT.Nf; column > 0; column--) {
        PHI_filt_window[row1][row2][column] = PHI_filt_window[row1][row2][column - 1];
      }
      PHI_filt_window[row1][row2][0] = PHI_filt[row1][row2];
    }
  }

  
   Print2DBuffer("u_filt",u_filt);
  for (int row = 0; row < FLAG.lz; row++) {
    for (int column = FILT.Nf; column > 0; column--) {
        u_filt_window[row][column] = u_filt_window[row][column - 1];
    }
    u_filt_window[row][0] = u_filt[row][0];
  }
  Print2DBuffer("PHI_filt",PHI_filt);
  Print3DBuffer("PHI_filt_window",PHI_filt_window);
  Print2DBuffer("u_filt",u_filt);
  Print2DBuffer("u_filt_window",u_filt_window);
}

void RCAC(int kk, float* u_in,float* z_in, float* yp_in, float* r_in){
  if(kk == 1){
    CalculateRegSize(); // Calculates the size of the gain matrix based on controller type
    InitializeBuffers(); // Initalized buffers to zero
    InitializeWindowBuffers();
  }
  
  UpdateBuffer(u_in, z_in, yp_in, r_in);
  std::vector<std::vector<float>> PHI = BuildRegressor();
  updateWindowBuffer(PHI);
  FilterSignals();
  // PrintBuffers();
  // Print3DBuffer("PHI_window", PHI_window);
  // Print2DBuffer("u_window", u_window);
  // Print2DBuffer("z_window", z_window);
  // Print2DBuffer("PHI",PHI);
  // Print3DBuffer("PHI_filt_window", PHI_filt_window);
  // Print2DBuffer("u_filt_window", u_filt_window);
  
}


void setup(){
    Serial.begin(9600);
    delay(100);
}
void loop(){
  
  u_in[0] = u_in[0] + 0.1;
  u_in[1] = u_in[1] + 0.1;

  z_in[0] = z_in[0] + 0.1;
  z_in[1] = z_in[1] + 0.1;

  yp_in[0] = yp_in[0] + 0.1;
  yp_in[1] = yp_in[1] + 0.1;

  r_in[0] = r_in[0] + 0.1;
  r_in[1] = r_in[1] + 0.1;

  RCAC(kk, u_in, z_in, yp_in, r_in);

  kk = kk + 1;
  delay(1000);
}

//----------------------------------------------------------------------------
//                    RCAC Functions -> Tested
//----------------------------------------------------------------------------



// Function to update the buffer with new values. Must take in an array
void UpdateBuffer(float* u_in, float* z_in, float* yp_in, float* r_in) {
    const int Nc = FLAG.Nc;
    const int lz = FLAG.lz;
    const int lu = FLAG.lu;
    const int ly = FLAG.ly;
    const int RegZ = FLAG.RegZ;

    // Update u_h
    for (int row = 0; row < lu; row++) {
        for (int column = Nc - 1; column > 0; column--) {
            u_h[row][column] = u_h[row][column - 1];
        }
        u_h[row][0] = (lu == 1) ? u_in[0] : u_in[row];
    }

    // Update z_h and r_h
    for (int row = 0; row < lz; row++) {
        for (int column = Nc; column > 0; column--) {
            z_h[row][column] = z_h[row][column - 1];
            r_h[row][column] = r_h[row][column - 1];
        }
        z_h[row][0] = (lz == 1) ? z_in[0] : z_in[row];
        r_h[row][0] = (lz == 1) ? r_in[0] : r_in[row];
    }

    // Update yp_h based on RegZ
    if (RegZ == 1) {
        for (int row = 0; row < lz; row++) {
            for (int column = Nc; column > 0; column--) {
                yp_h[row][column] = yp_h[row][column - 1];
            }
            yp_h[row][0] = (lz == 1) ? z_in[0] : z_in[row];
        }
    } else {
        for (int row = 0; row < ly; row++) {
            for (int column = Nc; column > 0; column--) {
                yp_h[row][column] = yp_h[row][column - 1];
            }
            yp_h[row][0] = (ly == 1) ? yp_in[0] : yp_in[row];
        }
    }
    for (int row = 0; row < lz; row ++){
      intg[row] = intg[row] + z_h[row][0];
    }
} // Tested and Works

// CalculateRegSize function
void CalculateRegSize() {
    if (FLAG.ControllerType == "PID") {
        ltheta = 2;
    } else if (FLAG.ControllerType == "PI") {
        ltheta = 1;
    } else if (FLAG.ControllerType == "P") {
        ltheta = 0;
    } else if (FLAG.ControllerType == "FSF") {
        ltheta = FLAG.lx;
    } else if (FLAG.ControllerType == "FSFI") {
        ltheta = (FLAG.lx+1 + FLAG.ly+1)-1;
    }
} // Tested and Works

// Function to initialize buffers and set up the system
void InitializeBuffers() {
    const int Nc = FLAG.Nc;
    const int lz = FLAG.lz;
    const int lu = FLAG.lu;
    const int ly = FLAG.ly;
    const int RegZ = FLAG.RegZ;

    // Initialize buffers
    u_h = new float*[lu];
    z_h = new float*[lz];
    r_h = new float*[lz];
    if (RegZ) {
        yp_h = new float*[lz];
    } else {
        yp_h = new float*[ly];
    }
    intg = new float[lz];

    // Allocate and initialize memory for each buffer
    for (int rows = 0; rows < lu; rows++) {
        u_h[rows] = new float[Nc];
        for (int cols = 0; cols < Nc; cols++) {
            u_h[rows][cols] = 0.0; // Initialize to zero
        }
    }
    for (int rows = 0; rows < lz; rows++) {
        z_h[rows] = new float[Nc + 1];
        r_h[rows] = new float[Nc + 1];
        for (int cols = 0; cols < Nc + 1; cols++) {
            z_h[rows][cols] = 0.0;
            r_h[rows][cols] = 0.0;
        }
    }

    // Allocate and initialize memory for yp_h based on RegZ
    if (RegZ) {
        for (int rows = 0; rows < lz; rows++) {
            yp_h[rows] = new float[Nc + 1];
            for (int cols = 0; cols < Nc + 1; cols++) {
                yp_h[rows][cols] = 0.0;
            }
        }
    } else {
        for (int rows = 0; rows < ly; rows++) {
            yp_h[rows] = new float[Nc + 1];
            for (int cols = 0; cols < Nc + 1; cols++) {
                yp_h[rows][cols] = 0.0;
            }
        }
    }

    for (int i = 0; i < lz; i++) {
        intg[i] = 0.0; // Initialize to zero
    }
} // Tested and Works

// Function to print buffers for debugging
void PrintBuffers() {
    Serial.println("Control Buffer (u_h):");
    for (int i = 0; i < FLAG.lu; i++) {
        Serial.print("u_h["); Serial.print(i); Serial.print("]: ");
        for (int j = 0; j < FLAG.Nc; j++) {
            Serial.print(u_h[i][j]); Serial.print(" ");
        }
        Serial.println();
    }

    Serial.println("Performance Buffer (z_h):");
    for (int i = 0; i < FLAG.lz; i++) {
        Serial.print("z_h["); Serial.print(i); Serial.print("]: ");
        for (int j = 0; j < FLAG.Nc + 1; j++) {
            Serial.print(z_h[i][j]); Serial.print(" ");
        }
        Serial.println();
    }

    Serial.println("Command Buffer (r_h):");
    for (int i = 0; i < FLAG.lz; i++) {
        Serial.print("r_h["); Serial.print(i); Serial.print("]: ");
        for (int j = 0; j < FLAG.Nc + 1; j++) {
            Serial.print(r_h[i][j]); Serial.print(" ");
        }
        Serial.println();
    }

    Serial.println("Regressor Buffer (yp_h):");
    for (int i = 0; i < (FLAG.RegZ ? FLAG.lz : FLAG.ly); i++) {
        Serial.print("yp_h["); Serial.print(i); Serial.print("]: ");
        for (int j = 0; j < FLAG.Nc + 1; j++) {
            Serial.print(yp_h[i][j]); Serial.print(" ");
        }
        Serial.println();
    }

    Serial.println("Integration Buffer (intg):");
    for (int i = 0; i < FLAG.lz; i++) {
        Serial.print("intg["); Serial.print(i); Serial.print("]: ");
        Serial.print(intg[i]);
        Serial.println();
    }
}

void InitializeWindowBuffers() {
    // Determine dimensions
    int nf_end = 5;
    int pc = nf_end;
    int pn = pc + FILT.Nf + FLAG.Nc;

    // Allocate dynamic buffers
    PHI_window.resize(FLAG.lu, std::vector<std::vector<float>>(ltheta, std::vector<float>(pc + pn, 0.0f)));
    u_window.resize(FLAG.lu, std::vector<float>(pn, 0.0f));
    z_window.resize(FLAG.lz, std::vector<float>(pn, 0.0f));

    PHI_filt_window.resize(FLAG.lz, std::vector<std::vector<float>>(ltheta, std::vector<float>(FILT.Nf, 0.0f)));
    u_filt_window.resize(FLAG.lz, std::vector<float>(FILT.Nf, 0.0f));
    theta_k.resize(ltheta, 0.0f);
    // Initialize Eigen matrix
    P_k = Eigen::MatrixXf::Identity(ltheta, ltheta) / FLAG.R0;
}


void Print3DBuffer(const char *name, const std::vector<std::vector<std::vector<float>>> &buffer) {
    Serial.println(name);
    for (size_t i = 0; i < buffer.size(); i++) {
        Serial.print("Layer ");
        Serial.println(i);
        for (size_t j = 0; j < buffer[i].size(); j++) {
            Serial.print("[");
            for (size_t k = 0; k < buffer[i][j].size(); k++) {
                Serial.print(buffer[i][j][k]);
                if (k < buffer[i][j].size() - 1) Serial.print(", ");
            }
            Serial.println("]");
        }
    }
    Serial.println();
}

void Print2DBuffer(const char *name, const std::vector<std::vector<float>> &buffer) {
    Serial.println(name);
    for (size_t i = 0; i < buffer.size(); i++) {
        Serial.print("[");
        for (size_t j = 0; j < buffer[i].size(); j++) {
            Serial.print(buffer[i][j]);
            if (j < buffer[i].size() - 1) Serial.print(", ");
        }
        Serial.println("]");
    }
    Serial.println();
}

// Function to print a 1D buffer
void Print1DBuffer(const char* name, const std::vector<float>& buffer) {
    Serial.println(name);
    Serial.print("[");
    for (size_t i = 0; i < buffer.size(); i++) {
        Serial.print(buffer[i]);
        if (i < buffer.size() - 1) Serial.print(", ");
    }
    Serial.println("]");
    Serial.println();
}
