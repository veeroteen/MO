//main.cpp
#include <iostream>
#include "LR2.h"
#include <iomanip>
#include <fstream>
#include "LR4.h"

int main()
{
   std::vector<double> C = { 2,4,2,6,2,3 };
   std::vector<double> a = { -3,-6,2,6,-3,8 };
   std::vector<double> b = { 6,-8,-8,8,-4,-1 };

   std::vector<std::pair<double, double>> area =
   {
      {-10,10},
      {-10,10}
   };

   auto fun = providerFUN(C, a, b, area);
   std::ofstream file("A1.csv");
   //std::vector<double> res = randomSearch(fun.first,1E-1,0.9,area,file);
   std::vector<double> res = A1(fun.first, fun.second, 1E-4, 1E-6, 100, area, file);
   file.close();
   std::cout << std::setprecision(8) << "(" << res[0] << "; " << res[1] << ")" << "\t" << -fun.first(res) << "\neval count: " << ccount << std::endl;

   
   return 0;
}


int old()
{
   auto a = providerRB();
   std::ofstream file("MSGRB.csv");

   std::vector<double> x = { 0,-1 };
   std::vector<double> res = MSGFR(a.first, a.second, x, 1E-5, 1E-7, 100000, file);
   std::cout << std::setprecision(16) << "(" << res[0] << "; " << res[1] << ")" << std::endl;

   return 0;
}