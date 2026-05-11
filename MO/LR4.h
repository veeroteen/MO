//LR4.h
#pragma once
#include <functional>
#include <utility>
#include <cmath>
#include <numbers>
#include <iostream>
#include <iomanip>
#include <array>
#include <algorithm>
#include <random>
#include "LR3.h"
#include <streambuf>

size_t ccount = 0;

struct Dot
{
   std::vector<double> cords;
   double value = 0;
   Dot(std::vector<double> &cords) : cords(cords){};
   bool operator >(Dot &b)
   {
      return value > b.value ? true : false;
   
   }
};



class pokeAtTheSky
{
   std::vector<std::uniform_real_distribution<>> random;
   std::random_device rd;
   std::mt19937 gen;
   std::function<double(std::vector<double> &)> fun;
public:
   pokeAtTheSky(std::vector<std::pair<double,double>> &area,std::function<double(std::vector<double>&)> fun) : fun(fun)
   {
      auto a = rd();
      std::cout << "seed: " << a << std::endl;
      gen = std::mt19937(a);
      for(auto a : area)
      {
         random.emplace_back(a.first, a.second);
      }
   }

   Dot poke()
   {
      std::vector<double> cords(random.size(),0);
      for(size_t i = 0; i < cords.size();i++)
      {
         cords[i] = random[i](gen);
      }
      Dot res(cords);
      return res;
   }
};

auto providerFUN(std::vector<double> &C, std::vector<double> &a, std::vector<double> &b, std::vector<std::pair<double,double>> &area)
{
      
   auto fun = std::pair
   {
      [&](const std::vector<double> &x)
      {
         double buff = 0;
         for (size_t i = 0; i < C.size(); i++)
         {
            buff += C[i] / (1 + std::pow(x[0] - a[i], 2) + std::pow(x[1] - b[i], 2));
         }
         ccount++;
         return -buff;
      },
      [&](const std::vector<double> &x)
      {
         double buff =
            -log(x[0] - area[0].first -1 ) 
            -log(-x[0] + area[0].second +1) 
            -log(x[1] - area[1].first -1) 
            -log(-x[1] + area[1].second +1);
         return isnan(buff) ? std::numeric_limits<double>::infinity() : 0;

      }
   };

   return fun;
}

std::vector<double> randomSearch
(
   std::function<double(const std::vector<double> &)> fun,
   double eps, double P, std::vector<std::pair<double, double>> &area, std::ostream &stream
)
{
   pokeAtTheSky radomizer(area, fun);
   double Veps = std::pow(eps, area.size());
   double V = 1;
   for(auto a : area)
   {
      V *= (a.second - a.first);
   }
   double Peps = Veps / V;
   size_t N = std::log(1 - P) / std::log(1 - Peps);
   std::cout << "dot count: " << N << std::endl;
   Dot res = radomizer.poke();
   res.value = fun(res.cords);
   stream << res.cords[0] << "\t" << res.cords[1] << std::endl;
   for(size_t i = 0; i < N; i++)
   {
      
      Dot buff = radomizer.poke();
      buff.value = fun(buff.cords);
      if(res > buff)
      {
         res = buff;
      }
      stream << res.cords[0] << "\t" << res.cords[1] << std::endl;
   }

   return(res.cords);
}

std::vector<double> A1
(
   std::function<double(const std::vector<double> &)> fun, std::function<double(const std::vector<double> &)> barf,
   double eps, double leps,size_t m, std::vector<std::pair<double, double>> &area, std::ostream &stream
)
{
   pokeAtTheSky radomizer(area, fun);
   Dot curD = radomizer.poke();
   curD.cords = BarrierMeth(fun, barf, curD.cords, eps, leps, 1000, stream);
   curD.value = fun(curD.cords);
   stream << std::setprecision(16) << curD.cords[0] << "\t" << curD.cords[1] << std::endl;
   size_t counter = 0;
   while(counter < m)
   {
      Dot nextD = radomizer.poke();
      nextD.cords = BarrierMeth(fun, barf, nextD.cords, eps, leps, 1000, nullStream());
      nextD.value = fun(nextD.cords);

      if(curD > nextD)
      {
         curD = nextD;
         stream<< std::setprecision(16) << nextD.cords[0] << "\t" << nextD.cords[1] << std::endl;
      }
      else 
      {
         counter++;
      }
   }
   return curD.cords;
}

std::vector<double> A2
(
   std::function<double(const std::vector<double> &)> fun, std::function<double(const std::vector<double> &)> barf,
   double eps, double leps, size_t m, std::vector<std::pair<double, double>> &area, std::ostream &stream
)
{
   pokeAtTheSky radomizer(area, fun);
   Dot curD = radomizer.poke();
   std::vector<double> cur = BarrierMeth(fun, barf, curD.cords, eps, leps, 1000, nullStream());
   curD.cords = cur;
   curD.value = fun(cur);
   stream << std::setprecision(16) << curD.cords[0] << "\t" << curD.cords[1] << std::endl;
   size_t counter = 0;
   while (counter < m)
   {

      Dot nextD = radomizer.poke();
      
      nextD.value = fun(nextD.cords);

      if (curD > nextD)
      {
         nextD.cords = BarrierMeth(fun, barf, nextD.cords, eps, leps, 1000, nullStream());
         nextD.value = fun(nextD.cords);
         curD = nextD;
         stream << std::setprecision(16) << nextD.cords[0] << "\t" << nextD.cords[1] << std::endl;
      }
      else
      {
         counter++;
      }
   }
   return curD.cords;
}

void normilize(std::vector<double> &dir)
{
   double buff = 0;
   for (auto &a : dir)
   {
      buff += a * a;
   }
   buff = sqrt(buff);
   for (auto &a : dir)
   {
      a /= buff;
   }
}

bool check(std::vector<std::pair<double, double>> &area, std::vector<double> dot)
{
   if(area[0].first > dot[0] || area[0].second < dot[0] || area[1].first > dot[1] || area[1].second < dot[1])
   {
      return false;
   }
   return true;
}


std::vector<double> A3
(
   std::function<double(const std::vector<double> &)> fun, std::function<double(const std::vector<double> &)> barf,
   double eps, double leps, size_t m, std::vector<std::pair<double, double>> &area, std::ostream &stream
)
{
   pokeAtTheSky radomizer(area, fun);
   Dot curD = radomizer.poke();
   std::vector<double> cur = BarrierMeth(fun, barf, curD.cords, eps, leps, 1000, nullStream());
   curD.cords = cur;
   curD.value = fun(cur);
   stream << std::setprecision(16) << curD.cords[0] << "\t" << curD.cords[1] << std::endl;
   size_t counter = 0;

   while (counter < m)
   {
      Dot dir = radomizer.poke();
      normilize(dir.cords);
      std::function<double(double)> layer = [&](double a)
         {
            std::vector<double> buff = { dir.cords[0] * a + curD.cords[0] ,dir.cords[1] * a + curD.cords[1] };
            double b = -fun(buff) + barf(buff);
            return b;
         };
      std::pair<double, double> interval = findInterval(layer, 10,5);
      double a = Gold(interval, leps, layer).second + 0.5;
      std::vector<double> dot = { dir.cords[0] * (a + 1) + curD.cords[0] ,dir.cords[1] * (a+1) + curD.cords[1] };
      Dot nextD(dot);
      if(!check(area,dot))
      {
         continue;
      }
      nextD.cords = BarrierMeth(fun, barf, nextD.cords, eps, leps, 1000, nullStream());
      nextD.value = fun(nextD.cords);
      if (curD > nextD)
      {
         curD = nextD;
         stream << std::setprecision(16) << nextD.cords[0] << "\t" << nextD.cords[1] << std::endl;
      }
      else
      {
         counter++;
      }
   }
   return curD.cords;
}