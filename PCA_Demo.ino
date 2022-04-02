#include <BasicLinearAlgebra.h>

using namespace BLA;

#define PCA_PC  (2)
#define PCA_OBS (34)
#define PCA_VAR (3)

void setup() {

  Serial.begin(9600);
  while (!Serial) {
  }

  float rawIn[] = { -17.07, -8.95, -81.16, -14.58, 11.67, -71.91, -6.25, 30.88, -61.28, 7.26, 46.89, -52.42, 24.9, 59.44, -44.36, 46.97, 66.09, -37.11, 68.7, 66.49, -33.79, 91.98, 61.72, -32.31, 113.77, 51.64, -33.7, 131.82, 35.93, -37.33, 145.41, 17.1, -44.45, 153.43, -6.18, -52.76, 156.9, -28.85, -62.21, 152.26, -51.37, -73.21, 142.07, -70.57, -83.02, 126.54, -85.68, -93.74, 106.84, -95.85, -101.66, 84.14, -100.05, -107.28, 60.31, -97.55, -110.5, 37.46, -89.7, -110.58, 16.28, -75.81, -107.48, -0.71, -58.14, -102.89, -11.02, -37.88, -93.26, -16.64, -16.05, -84.42, -15.16, 6.53, -73.46, -7.9, 27.6, -63.15, 5.75, 44.37, -52.95, 24.04, 57.18, -44.29, 45.2, 64.9, -38.15, 67.72, 67.29, -35.45, 90.41, 62.17, -32.69, 111.89, 52.59, -33.56, 130.63, 36.67, -38.04, 144.9, 17.15, -44.76, 153.26, -5.36, -51.79, 156.23, -28.04, -61.25};


  BLA::Matrix<PCA_OBS, PCA_VAR> datamatrix;
  CreateMatrix(datamatrix, rawIn);                // Fill the matrix
  //Serial << datamatrix.Submatrix<4,PCA_VAR>(0,0);

  BLA::Matrix<1, PCA_VAR> centroid_;
  Centroid(datamatrix, centroid_ );               // Get the centroid
  //Serial << centroid_;

  // Subtract the centroid
  for ( int i = 0; i < PCA_VAR; i++) {
    datamatrix.Submatrix<PCA_OBS, 1>(0, i) -=  centroid_(0, i);
  }

  BLA::Matrix<PCA_OBS, PCA_PC> T; // T: scores or transformed coordinates.    dim: Datapoints x PC
  BLA::Matrix<PCA_VAR, PCA_PC> P; // P: loadings or coefficients.             dim: Variable count x PC
  int iterations = Nipals(datamatrix, T, P);                           // Run NIPALS algorithm

  Serial.print("iterations: ");
  Serial.println(iterations);
  //  Serial.println("T scores");
  //  Serial << T << '\n';
  //  Serial.println("P Loadings");
  //  Serial << P << '\n';


  // How to convert new data into the PC-space. (assumes data is already centered).
  BLA::Matrix<PCA_OBS, PCA_VAR> newdata = datamatrix;
  for ( int i = 0; i < newdata.Rows; i++) {
    float angle = AngleInPCSpace(P, newdata.Submatrix<1, PCA_VAR>(i, 0));
    Serial.println(angle);
  }
}

void loop() {
  // put your main code here, to run repeatedly:

}

float AngleInPCSpace(BLA::Matrix<PCA_VAR, PCA_PC> & my_P, BLA::Matrix<1, PCA_VAR> instance) {
  float pc1_score = (~my_P.Submatrix<PCA_VAR, 1>(0, 0) * ~instance)(0, 0); // first col in P transposed * a row in newdata transposed;
  float pc2_score = (~my_P.Submatrix<PCA_VAR, 1>(0, 1) * ~instance)(0, 0); // second col in P transposed * a row in newdata transposed;
  return atan2(pc1_score, pc2_score);
}

void Centroid(BLA::Matrix<PCA_OBS, PCA_VAR> data, BLA::Matrix<1, PCA_VAR> & centroid ) {
  centroid.Fill(0);
  for (int i = 0; i < PCA_OBS; i++) {
    float weight = 1.0 / (i + 1);
    centroid += (data.Submatrix<1, PCA_VAR>(i, 0) - centroid) * weight;
  }
}

void CreateMatrix(class BLA::Matrix<PCA_OBS, PCA_VAR> & mat, float datas[]) {
  int _col = mat.Cols;
  int _row = mat.Rows;
  //BLA::Matrix<PCA_OBS, PCA_VAR> mat;
  for (int i = 0; i < _row; i++) {
    for (int j = 0; j < _col; j++) {
      mat(i, j) = datas[i * _col + j];
    }
  }
}

int Nipals(BLA::Matrix<PCA_OBS, PCA_VAR> Xh, BLA::Matrix<PCA_OBS, PCA_PC> & T, BLA::Matrix<PCA_VAR, PCA_PC> & P) {
  // The function assumes already centered data.
  int it_total = 0;
  const int it = 1000; // Maximum iterations per PC
  const float tol = 1e-4; // Tolerance

  //T.Fill(0);
  //P.Fill(0);
  BLA::Matrix<PCA_OBS, 1> thnew;

  int nr = 0;
  for (int h = 0; h < PCA_PC; h++) {  // Loop for each principal component
    BLA::Matrix<PCA_OBS, 1> th =  Xh.Column(0);
    BLA::Matrix<PCA_VAR, 1> ph;
    ph.Fill(0);
    bool ende = false;
    while (!ende) {   // loop until the change in p (loadings) gets small.
      nr += 1;      
      ph = ~Xh * th / (~th * th)(0);
      ph = ph / Norm(ph);      
      thnew = Xh * ph / (~ph * ph)(0);
      float prec = (~(thnew - th) * (thnew - th))(0);
      th = thnew;
      if (prec <= tol * tol) {
        ende = true;
      }
      else if (it <= nr) {
        ende = true;
        Serial.println("PCA: no convergence");
      }
    }
    Xh = Xh - th * ~ph;
    T.Submatrix<PCA_OBS, 1>(0, h) = th;
    P.Submatrix<PCA_VAR, 1>(0, h) = ph;

    it_total += nr;
    nr = 0;
  }
  return it_total;
}


void printMatrix(BLA::Matrix<20, 1> mat) {
  for (int i = 0; i < mat.Rows; ++i)
  {
    for (int j = 0; j < mat.Cols; ++j)
    {
      Serial.print(mat(i, j));
      Serial.print(" , ");
    }
  }
}

float Norm(BLA::Matrix<PCA_VAR> sample) {
  float sum = 0;
  for (int i = 0; i < PCA_VAR; i++) {
    sum += sample(i) * sample(i);
  }
  return sqrt(sum);
}
