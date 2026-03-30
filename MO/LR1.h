#pragma once
#include <functional>
#include <utility>
#include <cmath>
#include <numbers>
#include <iostream>
#include <iomanip>
size_t counter = 0;
//Дихотомия
std::pair<double, double> Dih(std::pair<double, double> interval, double epsilon, std::function<double(double)> f, std::ostream &stream)
{
   size_t n = std::log2((interval.second - interval.first) / epsilon) + 1 ;
   double delta = epsilon / 2;
   std::cout << n << std::endl;
   std::pair<double, double> newIntr;
   for (size_t i = 0; i < n; i++)
   {
      newIntr=
      {
            (interval.first + interval.second - delta) / 2,
            (interval.first + interval.second + delta) / 2
      };
      stream << newIntr.first << "\t" << newIntr.second << "\t" << f(newIntr.first) << "\t" << f(newIntr.second) << "\t" << interval.first << "\t" << interval.second << "\n";
      if (f(newIntr.first) < f(newIntr.second))
      {
         interval.second = newIntr.first;
      }
      else if ((f(newIntr.first) > f(newIntr.second)))
      {
         interval.first = newIntr.second;
      }
      else
      {
         interval = newIntr;
      }

   }
   stream << newIntr.first << "\t" << newIntr.second << "\t" << f(newIntr.first) << "\t" << f(newIntr.second) << "\t" << interval.first << "\t" << interval.second << "\n";
   return interval;
}
//Золотое сечение
std::pair<double, double> Gold(std::pair<double, double> interval, double epsilon, std::function<double(double)> f)
{
   
   constexpr double phi = std::numbers::phi_v<long double>;
   size_t n = std::log((interval.second - interval.first) / epsilon) / std::log(phi) + 1;
   double delta = interval.second - interval.first;
   std::pair<double, double> newIntr
   {
      interval.second - delta/ phi,
      interval.first + delta/ phi
   };
   std::pair<double, double> funV
   {
      f(newIntr.first),
      f(newIntr.second)
   };
   std::cout << n << std::endl;
   for (size_t i = 0; i < n; i++)
   {
      if(funV.first < funV.second)
      {
         interval.second = newIntr.second;
         delta = interval.second - interval.first;

         newIntr.second = newIntr.first;
         newIntr.first = interval.second - delta / phi;
         funV.second = funV.first;
         funV.first = f(newIntr.first);
      }
      else
      {
         interval.first = newIntr.first;
         delta = interval.second - interval.first;

         newIntr.first = newIntr.second;
         newIntr.second = interval.first + delta / phi;
         funV.first = funV.second;
         funV.second = f(newIntr.second);
      }
   }
   return interval;
}
