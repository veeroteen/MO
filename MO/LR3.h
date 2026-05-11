//LR3.h
#pragma once
#include <functional>
#include <utility>
#include <cmath>
#include <numbers>
#include <iostream>
#include <iomanip>
#include <array>
#include <algorithm>
#include "LR2.h"
size_t count = 0;
auto providerPA()
{
   auto a = std::pair
   (
      [](const std::vector<double> &x)
      {
         count++;
         return 10*pow(x[1] - x[0], 2) + pow(x[1] ,2);
      },
      [](const std::vector<double> &x)
      {
         return std::max(0.0,pow(1 - x[0] - x[1],2));
      }
   );
   return a;
}

auto providerPB()
{
   auto a = std::pair
   (
      [](const std::vector<double> &x)
      {
         count++;
         return 10 * pow(x[1] - x[0], 2) + pow(x[1], 2);
      },

      [](const std::vector<double> &x)
      {
         return abs(pow(x[0] + x[1] - 2,2));
      }
   );
   return a;
}

std::vector<double> PenaltyMeth
(
   std::function<double(const std::vector<double> &)> fun,
   std::function<double(const std::vector<double> &)> penf,
   const std::vector<double> &x0,
   double eps, double leps, size_t maxI, std::ostream &stream
)
{
   
   double r = 5,q = 10;
   std::vector<double> x = x0;
   std::vector<double> xp = x;
   auto layer = [&](const std::vector<double> &cur)
      {
         return fun(cur) + r * penf(cur);
      };
   double v = 0;
   size_t i = 0;
   stream << x[0] << "\t" << x[1] << std::endl;
   for(; i < maxI; i++)
   {
      x = Gauss(layer, x, eps, leps, maxI, std::cout);
      v = penf(x);
      r *= q;
      if (v < eps && mod(x, xp) < eps)
      {
         break;
      }
      stream << x[0] << "\t" << x[1] << std::endl;
      xp = x;
   }
   std::cout << "iteration count: " << i+1 << std::endl;
   std::cout << "evaluation count: " << count << std::endl;
   return x;
}

auto providerBAl()
{
   auto a = std::pair
   (
      [](const std::vector<double> &x)
      {
         count++;
         return 10 * pow(x[1] - x[0], 2) + pow(x[1], 2);
      },
      [](const std::vector<double> &x)
      {
         return -log(-1 + x[0] + x[1]);
      }
   );
   return a;
}

auto providerBAd()
{
   auto a = std::pair
   (
      [](const std::vector<double> &x)
      {
         count++;
         return 10 * pow(x[1] - x[0], 2) + pow(x[1], 2);
      },
      [](const std::vector<double> &x)
      {
         auto a = 1 / (-1 + x[0] + x[1]);
         return a <= 0 ? std::numeric_limits<double>::infinity() : a;
      }
   );
   return a;
}

std::vector<double> BarrierMeth
(
   std::function<double(const std::vector<double> &)> fun,
   std::function<double(const std::vector<double> &)> barf,
   const std::vector<double> &x0,
   double eps, double leps, size_t maxI, std::ostream &stream
)
{
   double u = 1, q = 0.05;
   std::vector<double> x = x0;
   std::vector<double> xp = x;
   auto layer = [&](const std::vector<double> &cur)
      {
         return fun(cur) + u * barf(cur);
      };
   size_t i = 0;
   stream << x[0] << "\t" << x[1] << std::endl;
   for (; i < maxI; i++)
   {
      x = Gauss(layer, x, eps, leps, maxI, std::cout);
      
      if (u < eps && mod(x, xp) < eps)
      {
         break;
      }
      u *= q;
      stream << x[0] << "\t" << x[1] << std::endl;
      xp = x;
   }
   return x;
}