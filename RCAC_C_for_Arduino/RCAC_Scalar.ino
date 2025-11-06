#include "arduino_secrets.h"

struct RCACFLAG {
  static const bool RegZ = true;
  static const float R0 = 10000.0;
  static const float Rz = 1.0;
  static const bool FeedForward = false;
  static const float lambda = 1.0;
  static const int Nc = 4;
  static const bool Integrator = false;
  String ControllerType = "PID";
  static const bool FIR = false;
  static const int window = 15;
} FLAG;

struct RCACFILT {
  static const int Nf = 5;
  static const float Nu = 1;
} FILT;

static const int lu = 1;
static const int lz = 1;
static const int ly = 1;
static const int lx = 1;
static int ltheta;

// H buffers
static float* u_h;
static float* z_h;
static float* r_h;
static float* yp_h;

// Window buffers
static float** PHI_window;
static float* u_window;
static float* z_window;

// Filt window buffers
static float** PHI_filt_window;
static float* u_filt_window;

static float intg = 0.0;

static float* PHI;
// A constant 3x3 matrix P_k
static float P_k[3][3] = {
  { 1.0 / FLAG.R0, 0.0, 0.0 },
  { 0.0, 1.0 / FLAG.R0, 0.0 },
  { 0.0, 0.0, 1.0 / FLAG.R0 }
};

void RCAC_Scalar(int k, float u_in, float z_in, float yp_in, float r_in) {
  if (k == 1) {
    ltheta = CalculateRegSize();
    InitializeWindowBuffers();
    InitalizeHBuffer();
  }
  UpdateHBuffers(u_in, z_in, yp_in, r_in);
  BuildRegressor(k);
  // updateWindowBuffer(u_in, z_in);
}

int CalculateRegSize() {
  int size;
  if (FLAG.ControllerType == "PID") {
    size = 3;
  } else if (FLAG.ControllerType == "PI") {
    size = 2;
  } else if (FLAG.ControllerType == "P") {
    size = 1;
  } else if (FLAG.ControllerType == "FSF") {
    size = 1;
  } else if (FLAG.ControllerType == "FSFI") {
    size = 2;
  }
  return size;
}

void InitalizeHBuffer() {
  // Note: u_h is allocated with size Nc-1 while others are size Nc.
  u_h = new float[FLAG.Nc - 1];
  for (int cols = 0; cols < FLAG.Nc - 1; cols++) {
    u_h[cols] = 0.0;
  }
  z_h = new float[FLAG.Nc];
  r_h = new float[FLAG.Nc];
  yp_h = new float[FLAG.Nc];
  for (int cols = 0; cols < FLAG.Nc; cols++) {
    z_h[cols] = 0.0;
    r_h[cols] = 0.0;
    yp_h[cols] = 0.0;
  }
}

void InitializeWindowBuffers() {
  // Determine dimensions for the window buffers
  static const int nf_end = 5;
  static const int pc = nf_end;
  static const int pn = pc + FILT.Nf + FLAG.Nc;

  // Allocate and initialize PHI_window (dimensions: (ltheta-1) x (pn-1))
  PHI_window = new float*[ltheta];
  for (int rows = 0; rows < ltheta - 1; rows++) {
    PHI_window[rows] = new float[pn];
    for (int cols = 0; cols < pn - 1; cols++) {
      PHI_window[rows][cols] = 0.0;
    }
  }

  // Allocate and initialize u_window  and z_window (size: pn-1)
  u_window = new float[pn];
  z_window = new float[pn];
  for (int cols = 0; cols < pn - 1; cols++) {
    u_window[cols] = 0.0;
    z_window[cols] = 0.0;
  }

  // Allocate and initialize PHI_filt_window (dimensions: (ltheta-1) x (2*ltheta-1))
  PHI_filt_window = new float*[ltheta];
  for (int rows = 0; rows < ltheta - 1; rows++) {
    PHI_filt_window[rows] = new float[2 * ltheta - 1];
    for (int cols = 0; cols < 2 * ltheta - 1; cols++) {
      PHI_filt_window[rows][cols] = 0.0;
    }
  }

  // Allocate and initialize u_filt_window (size: 2*ltheta-1)
  u_filt_window = new float[2 * ltheta];
  for (int cols = 0; cols < 2 * ltheta - 1; cols++) {
    u_filt_window[cols] = 0.0;
  }
}

// UpdateHBuffers: Shifts all history buffers and updates them with new inputs
void UpdateHBuffers(float u_in, float z_in, float yp_in, float r_in) {
  // Update u_h (size: FLAG.Nc - 1)
  int len_u = FLAG.Nc - 1;
  // Shift right: for i from (len_u-1) downto 1, assign u_h[i] = u_h[i-1]
  for (int i = len_u - 1; i > 0; i--) {
    u_h[i] = u_h[i - 1];
  }
  // Set the first element to new control input
  u_h[0] = u_in;

  // Update z_h (size: FLAG.Nc)
  int len_z = FLAG.Nc;
  for (int i = len_z - 1; i > 0; i--) {
    z_h[i] = z_h[i - 1];
  }
  z_h[0] = z_in;

  // Update r_h (size: FLAG.Nc)
  int len_r = FLAG.Nc;
  for (int i = len_r - 1; i > 0; i--) {
    r_h[i] = r_h[i - 1];
  }
  r_h[0] = r_in;

  // Update yp_h (size: FLAG.Nc)
  int len_yp = FLAG.Nc;
  for (int i = len_yp - 1; i > 0; i--) {
    yp_h[i] = yp_h[i - 1];
  }
  // Choose value based on FLAG.RegZ
  if (FLAG.RegZ == true) {
    yp_h[0] = z_in;
  } else {
    yp_h[0] = yp_in;
  }

  // Update the integrator variable with the new z_h value
  intg = intg + z_h[0];
}

void BuildRegressor(int k) {
  if (k == 1) {
    PHI = new float[ltheta];
  }
  if (FLAG.ControllerType == "PID") {

    PHI[0] = yp_h[0];
    PHI[1] = intg;
    PHI[2] = yp_h[0] - yp_h[1];
  } else if (FLAG.ControllerType == "PI") {
    PHI[0] = yp_h[0];
    PHI[1] = intg;
  } else if (FLAG.ControllerType == "P") {
    PHI[0] = yp_h[0];
  }
}

void updateWindowBuffer(float u_in, float z_in) {
  static const int nf_end = 5;
  static const int nf = FILT.Nf;
  static const int length_window = nf + nf_end;

  for (int column = length_window; column > 0; column--) {
    u_window[column] = u_window[column - 1];
  }
  u_window[0] = u_in;

  for (int column = length_window; column > 0; column--) {
    z_window[column] = z_window[column - 1];
  }
  z_window[0] = z_in;


  for (int row = 0; row < ltheta; row++) {
    for (int column = length_window; column > 0; column--) {
      PHI_window[row][column] = PHI_window[row][column - 1];
    }
    PHI_window[row][0] = PHI[row];
  }
}


void PrintBuffers() {
  // Print H Buffers
  Serial.println("=== H Buffers ===");

  Serial.print("u_h: ");
  for (int i = 0; i < FLAG.Nc - 1; i++) {
    Serial.print(u_h[i]);
    if (i < FLAG.Nc - 2)
      Serial.print(", ");
  }
  Serial.println();

  Serial.print("z_h: ");
  for (int i = 0; i < FLAG.Nc; i++) {
    Serial.print(z_h[i]);
    if (i < FLAG.Nc - 1)
      Serial.print(", ");
  }
  Serial.println();

  Serial.print("r_h: ");
  for (int i = 0; i < FLAG.Nc; i++) {
    Serial.print(r_h[i]);
    if (i < FLAG.Nc - 1)
      Serial.print(", ");
  }
  Serial.println();

  Serial.print("yp_h: ");
  for (int i = 0; i < FLAG.Nc; i++) {
    Serial.print(yp_h[i]);
    if (i < FLAG.Nc - 1)
      Serial.print(", ");
  }
  Serial.println();

  // Recalculate window dimensions for printing purposes
  static const int nf_end = 5;
  static const int pc = nf_end;
  static const int pn = pc + FILT.Nf + FLAG.Nc;

  // Print Window Buffers
  Serial.println("=== Window Buffers ===");
  Serial.println("PHI_window:");
  for (int i = 0; i < ltheta - 1; i++) {
    for (int j = 0; j < pn - 1; j++) {
      Serial.print(PHI_window[i][j]);
      if (j < pn - 2)
        Serial.print(", ");
    }
    Serial.println();
  }

  Serial.print("u_window: ");
  for (int i = 0; i < pn - 1; i++) {
    Serial.print(u_window[i]);
    if (i < pn - 2)
      Serial.print(", ");
  }
  Serial.println();

  Serial.print("z_window: ");
  for (int i = 0; i < pc - 1; i++) {
    Serial.print(z_window[i]);
    if (i < pc - 2)
      Serial.print(", ");
  }
  Serial.println();

  // Print Filt Window Buffers
  Serial.println("=== Filt Window Buffers ===");
  Serial.println("PHI_filt_window:");
  for (int i = 0; i < ltheta - 1; i++) {
    for (int j = 0; j < 2 * ltheta - 1; j++) {
      Serial.print(PHI_filt_window[i][j]);
      if (j < 2 * ltheta - 2)
        Serial.print(", ");
    }
    Serial.println();
  }

  Serial.print("u_filt_window: ");
  for (int i = 0; i < 2 * ltheta - 1; i++) {
    Serial.print(u_filt_window[i]);
    if (i < 2 * ltheta - 2)
      Serial.print(", ");
  }
  Serial.println();

  // // Print the constant matrix P_k
  // Serial.println("=== Constant Matrix P_k ===");
  // for (int i = 0; i < 3; i++) {
  //   for (int j = 0; j < 3; j++) {
  //     Serial.print(P_k[i][j], 6);
  //     if (j < 2)
  //       Serial.print(", ");
  //   }
  //   Serial.println();
  // }

  // Print the regressor vector PHI
  Serial.println("=== Regressor PHI ===");
  if (PHI != NULL) {
    for (int i = 0; i < ltheta; i++) {
      Serial.print(PHI[i]);
      if (i < ltheta - 1) {
        Serial.print(", ");
      }
    }
    Serial.println();
  } else {
    Serial.println("PHI is not allocated.");
  }
}

float k = 1;
float u = 1;
float z = 2;
float yp = 3;
float r = 4;

void setup() {
  Serial.begin(9600);
  delay(50);
  RCAC_Scalar(k, u, z, yp, r);
  PrintBuffers();
}

void loop() {
  // tested agianst matlab, matches results
  u = u + 0.1;
  z = z + 0.11;
  yp = yp + 0.12;
  r = r + 0.13;
  k = k + 1;
  RCAC_Scalar(k, u, z, yp, r);
  PrintBuffers();
  delay(1000);
}
