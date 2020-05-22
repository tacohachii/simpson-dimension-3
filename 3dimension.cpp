#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

const double mass = 1.0; //kg
const double rho_0 = 1.0; //kg.m-3
const double R = 1.0; //m
const double G = 6.67430e-11; //m+3.kg-1.s-2
const double PI = acos(-1.0); //3.14159265...

// n数をもっと増やしたいが，計算時間増大する
const int n_x = 100; //x軸の刻み数
const int n_y = 100; //y軸の刻み数
const int n_z = 100; //z軸の刻み数

const double x_min = 0.0 , x_max = PI ; //xの積分範囲
const double y_min = 0.0 , y_max = 1.0 ; //yの積分範囲
const double z_min = -1.0 , z_max = 1.0 ; //zの積分範囲

//被積分関数
double F(double x, double y, double z) 
{
  return y*sin(x)+z*sin(x); // ここを書き換える
}

// シンプソン法
double simpe3(const int n_x, const int n_y, const int n_z)
{
  int i, j, k;
  double h1, h2, h3, r, theta, phi;
  double v;

  double ***f = new double**[n_x+1]; // n_x+1 個分の領域を動的確保
  for (i= 0; i<=n_x; ++i) {
    f[i] = new double*[n_y+1]; // n_y+1 個分の領域を動的確保
    for (j= 0; j<=n_y; ++j) {
      f[i][j] = new double[n_z+1]; // n_z+1 個分の領域を動的確保
    }
  }

  for(i = 0; i <= n_x; i++){
    for (j = 0; j <= n_y; j++)  {
      for (k = 0; k <= n_z; k++)  {
        r = x_min + (x_max-x_min)/double(n_x) * double(i);
        theta = y_min + (y_max-y_min)/double(n_y) * double(j);
        phi = z_min + (z_max-z_min)/double(n_z) * double(k);
        f[i][j][k] = F(r,theta,phi);
      }
    }
  }
  h1 = (x_max-x_min) / double(n_x);
  h2 = (y_max-y_min) / double(n_y);
  h3 = (z_max-z_min) / double(n_z);

  // ---------------------　計算　---------------------

  //n_x+1 個分の領域を動的確保
  double *temp = new double[n_x+1];  

  //n_x+1 x n_y+1 個分の領域を動的確保
  double **temp2 = new double*[n_x+1];  
  for (i= 0; i<=n_x; ++i) temp2[i]= new double[n_y+1];

  for(i = 0; i <= n_x; i++){
    for(j = 0; j <= n_y; j++){
      v = - f[i][j][0] + f[i][j][n_z];
      for(k = 0 ; k < n_z - 1; k += 2)
        v += (2 * f[i][j][k] + 4 * f[i][j][k + 1]);
      temp2[i][j] = v;
    }
    v = - temp2[i][0] + temp2[i][n_y];
    for(j = 0; j < n_y - 1; j += 2) 
      v += (2 * temp2[i][j] + 4 * temp2[i][j + 1]);
    temp[i] = v;
  }
  v = - temp[0] + temp[n_x];
  for(i = 0; i < n_x - 1; i += 2) 
    v += (2 * temp[i] + 4 * temp[i + 1]);

  // 動的に確保した領域をそれぞれ解放
  for (i= 0; i<=n_x; ++i) delete[] temp2[i]; 
  delete [] temp2;
  delete [] temp;

  // ---------------------　計算　---------------------

  // 動的に確保した領域をそれぞれ解放  
  for (i= 0; i<=n_x; ++i) {
    for (j= 0; j<=n_y; ++j) {
      delete[] f[i][j];
    }
    delete[] f[i];
  }
  delete[] f;

  return v * h1 * h2 * h3 / 27.0;
}

int main(){
  // double ans = -1*mass*G*rho_0*simpe3(n_x, n_y, n_z); 
  double ans = simpe3(n_x, n_y, n_z); // 2 https://jp.mathworks.com/help/matlab/ref/integral3.html
  cout << std::fixed << std::setprecision(15) << ans << endl;
  return 0;
}

