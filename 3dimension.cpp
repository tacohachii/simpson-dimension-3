#include <iostream>
#include <cmath>
using namespace std;

const double mass = 1.0; //kg
const double rho_0 = 1.0; //kg.m-3
const double R = 1.0; //m
const double G = 6.67430e-11; //m+3.kg-1.s-2
const double PI = acos(-1.0); //3.14159265...

// n数をもっと増やしたいが，計算時間増大する
const int n_r = 100; //r軸の刻み数
const int n_theta = 100; //theta軸の刻み数
const int n_phi = 100; //phi軸の刻み数

const double r_0 = 0.0 , r_max = R ; //rの積分範囲
const double theta_0 = 0.0 , theta_max = PI ; //thetaの積分範囲
const double phi_0 = 0.0 , phi_max = 2*PI ; //phiの積分範囲

//被積分関数
double F(double r, double theta, double phi) 
{
  return r*sin(theta)*sin(phi);
}

// シンプソン法
double simpe3(double ***f, const int n_r, const int n_theta, const int n_phi, double h1, double h2, double h3)
{
  int i, j, k;
  double v;

  //n_r+1 個分の領域を動的確保
  double *temp = new double[n_r+1];  

  //n_r+1 x n_theta+1 個分の領域を動的確保
  double **temp2 = new double*[n_r+1];  
  for (i= 0; i<=n_r; ++i) temp2[i]= new double[n_theta+1];

  for(i = 0; i <= n_r; i++){
    for(j = 0; j <= n_theta; j++){
      v = - f[i][j][0] + f[i][j][n_phi];
      for(k = 0 ; k < n_phi - 1; k += 2)
        v += (2 * f[i][j][k] + 4 * f[i][j][k + 1]);
      temp2[i][j] = v;
    }
    v = - temp2[i][0] + temp2[i][n_theta];
    for(j = 0; j < n_theta - 1; j += 2) 
      v += (2 * temp2[i][j] + 4 * temp2[i][j + 1]);
    temp[i] = v;
  }
  v = - temp[0] + temp[n_r];
  for(i = 0; i < n_r - 1; i += 2) 
    v += (2 * temp[i] + 4 * temp[i + 1]);

  // 動的に確保した領域をそれぞれ解放
  for (i= 0; i<=n_r; ++i) delete[] temp2[i]; 
  delete [] temp2;
  delete [] temp;

  return v * h1 * h2 * h3 / 27.0;
}

int main(){
  double h1, h2, h3, r, theta, phi;
  int i, j, k;

  double ***f = new double**[n_r+1]; // n_r+1 個分の領域を動的確保
  for (i= 0; i<=n_r; ++i) {
    f[i] = new double*[n_theta+1]; // n_theta+1 個分の領域を動的確保
    for (j= 0; j<=n_theta; ++j) {
      f[i][j] = new double[n_phi+1]; // n_phi+1 個分の領域を動的確保
    }
  }

  for(i = 0; i <= n_r; i++){
    for (j = 0; j <= n_theta; j++)  {
      for (k = 0; k <= n_phi; k++)  {
        r = r_0 + (r_max-r_0)/double(n_r) * double(i);
        theta = theta_0 + (theta_max-theta_0)/double(n_theta) * double(j);
        phi = phi_0 + (phi_max-phi_0)/double(n_phi) * double(k);
        f[i][j][k] = F(r,theta,phi);
      }
    }
  }
  h1 = (r_max-r_0) / double(n_r);
  h2 = (theta_max-theta_0) / double(n_theta);
  h3 = (phi_max-phi_0) / double(n_phi);
  
  double ans;
  ans = -1 * mass * G * rho_0 * simpe3(f, n_r, n_theta, n_phi, h1, h2, h3);

  // 動的に確保した領域をそれぞれ解放  
  for (i= 0; i<=n_r; ++i) {
    for (j= 0; j<=n_theta; ++j) {
      delete[] f[i][j];
    }
    delete[] f[i];
  }
  delete[] f;

  cout << ans << endl;
  return 0;
}

