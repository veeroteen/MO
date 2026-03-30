#include <iostream>
#include "LR2.h"
#include <iomanip>
#include <fstream>

size_t cc = 0;
auto providerRB()
{
   auto a = std::pair
   (
      [](const std::vector<double> &x)
      {

         return 100 * pow((x[1] - x[0] * x[0]), 2) + pow((1 - x[0]), 2);
      },
      [](const std::vector<double> &x)
      {
         std::vector<double> grad(x.size(), 0);
         grad[0] = -400 * (x[1] - x[0] * x[0]) * x[0] - 2 * (1 - x[0]);
         grad[1] = 200 * (x[1] - x[0] * x[0]);
         return grad;
      }
   );
   return a;
}

auto providerBase()
{
   auto a = std::pair
   (
      [](const std::vector<double> &x)
      {

         return 100 * (x[1] - x[0]) * (x[1] - x[0]) + (1 - x[0]) * (1 - x[0]);
      },
      [](const std::vector<double> &x)
      {
         std::vector<double> grad(x.size(), 0);
         grad[0] = -200 * (x[1] - x[0]) - 2 * (1 - x[0]);
         grad[1] = 200 * (x[1] - x[0]);
         return grad;
      }
   );
   return a;
}

auto FUN()
{
   auto a = std::pair
   (
      [](const std::vector<double> &x)
      {

         return -(1 / (1 + pow((x[0] - 2) / 3, 2) + pow((x[1] - 2) / 3, 2)) + 3 / (1 + pow((x[0] - 1) / 1, 2) + pow((x[1] - 1) / 2, 2)));
      },
      [](const std::vector<double> &x)
      {
         std::vector<double> grad(x.size(), 0);
         grad[0] = -( - 1 * 2 * ((x[0] - 2) / 3) / pow((1 + pow((x[0] - 2) / 3, 2) + pow((x[1] - 2) / 3, 2)), 2) - 3 * 2 * (x[0] - 1) / pow((1 + pow((x[0] - 2) / 3, 2) + pow((x[1] - 2) / 3, 2)), 2));
         grad[1] = -( - 1 * 2 * ((x[1] - 2) / 3) / pow((1 + pow((x[0] - 2) / 3, 2) + pow((x[1] - 2) / 3, 2)), 2) - 3 * 2 * ((x[1] - 1) / 2) / pow((1 + pow((x[0] - 2) / 3, 2) + pow((x[1] - 2) / 3, 2)), 2));
         return grad;
      }
   );
   return a;
}

int main()
{
   auto a = FUN();
   std::ofstream file("Broiden.csv");

   counter = 0;
   std::vector<double> x = { 0,0 };
   std::vector<double> res = MSGFR(a.first, a.second, x, 1E-5, 1E-7, 100000, file);
   std::cout << std::setprecision(16) << "(" << res[0] << "; " << res[1] << ")" << std::endl;
   
}

