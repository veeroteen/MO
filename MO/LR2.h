#pragma once
#include <functional>
#include <utility>
#include <cmath>
#include <numbers>
#include <iostream>
#include <iomanip>
#include <array>
#include "LR1.h"
#include <cmath>


double scalar(const std::vector<double> &a, const std::vector<double> &b)
{
   double buff = 0;
   for (size_t i = 0; i < a.size(); i++)
   {
      buff += a[i] * b[i];
   }
   return buff;
}

double mod(const std::vector<double> &a, const std::vector<double> &b)
{
   std::vector<double> tmp(a);
   for (size_t i = 0; i < tmp.size(); i++)
   {
      tmp[i] -= b[i];
   }
   return sqrt(scalar(tmp, tmp));

}

std::pair<double, double> findInterval(std::function<double(double)> &fun, double h)
{
   double a0 = 0.0;
   double a1 = h;

   double f0 = fun(a0);
   double f1 = fun(a1);

   while (true)
   {
      double a2 = 2 * a1;
      double f2 = fun(a2);

      if (f2 >= f1)
         return std::pair<double, double>( a0, a2 );

      a0 = a1;
      f0 = f1;

      a1 = a2;
      f1 = f2;
   }
}

void setGrad(std::function<double(const std::vector<double> &)> fun, const std::vector<double> &x,std::vector<double> &grad)
{
   size_t n = x.size();
   grad.resize(n);
   double delta = 1E-8;
   for (size_t i = 0; i < n; i++)
   {
      std::vector<double> xp1 = x;
      std::vector<double> xp2 = x;
      std::vector<double> xm1 = x;
      std::vector<double> xm2 = x;

      xp1[i] += delta;
      xp2[i] += 2 * delta;
      xm1[i] -= delta;
      xm2[i] -= 2 * delta;

      grad[i] = (-fun(xp2) + 8 * fun(xp1) - 8 * fun(xm1) + fun(xm2)) / (12 * delta);
   }
}

std::pair<double,double> squareInt(std::pair<double, double> interval, double eps, std::function<double(double)> f)
{
   double a1, a2, a3;
   a1 = interval.first;
   a3 = interval.second;
   a2 = (a1 + a3) / 2;
   double f1, f2, f3;
   f1 = f(a1);
   f2 = f(a2);
   f3 = f(a3);
   double a = 0;

   while(fabs(a3-a1) > eps && f1 != f3 && fabs(f3 - f1) > eps)
   {
      
      double denom = ((a2 - a1) * (f2 - f3) - (a2 - a3) * (f2 - f1));
      double numer = (pow((a2 - a1),2)* (f2 - f3) - pow((a2 - a3),2) * (f2 - f1));
      
      
      double tmp = a2 - numer / (denom * 2);
      if(fabs(tmp-a) < eps)
      {
         break;
      }
      a = tmp;
      double fr = f(a);

      if(a < a2)
      {
         if(fr < f2)
         {
            a3 = a2;
            f3 = f2;
            a2 = a;
            f2 = fr;
         }
         else
         {
            a1 = a;
            f1 = fr;
         }
      }
      else
      {
         if(fr < f2)
         {
            a1 = a2;
            f1 = f2;
            a2 = a;
            f2 = fr;
         }
         else
         {
            a3 = a;
            f3 = fr;
         }
      }
      counter++;
   }
   return std::pair <double,double> (a2, a2);

}

std::vector<double> MSGFR
(
   std::function<double(const std::vector<double>&)> fun, 
   std::function<std::vector<double>(const std::vector<double>&)> gf,
   const std::vector<double> &x0,
   double eps,double leps,size_t maxI, std::ostream &stream
)
{
   std::vector<double> grad(x0.size());
   size_t k = 0;
   std::vector<double> direct(x0.size());
   std::vector<double> cur(x0);
   grad = gf(x0);
   for(size_t i = 0; i < direct.size();i++)
   {
      direct[i] = -grad[i];
   }
   std::function<double(double)> layer = [&fun,&cur,&direct](double lambda)
      {
         std::vector<double> arg(cur);
         for(size_t i = 0; i < arg.size();i++)
         {
            arg[i] += lambda * direct[i];
         }
         return fun(arg);
      };
   double funt = fun(x0);
   double gg = scalar(grad,grad);
   double a = 0;
   //stream << "x" << "\t" << "y" << "\t" << "f" << "\t" << "dirX" << "\t" << "dirY" << "\t" << "lambda" << "\t" << "dx" << "\t" << "dy" << "\t" << "df" << "\t" << "phi" << "\t" << "gradX" << "\t" << "gradY" << std::endl;

   do
   {
      for (auto &i : cur)
      {
         stream << std::setprecision(16) << i << "\t";
      }

      std::pair<double, double> interval = findInterval(layer, 0.5);
      a = squareInt(interval, leps, layer).second;
      if (a > 1.0) a = 1.0;
      if (a < 1e-8) a = 1e-8;

      double func = 0;
      std::vector<double> tmp(cur);
      while (true)
      {
         tmp = cur;
         for (size_t i = 0; i < cur.size(); i++)
         {
            tmp[i] += a * direct[i];
         }
         func = fun(tmp);
         if (func > funt)
         {
            a *= 0.8;
         }
         else
         {
            break;
         }
         if (a < 1e-12)
         {
            break;
         }
      }
      stream << funt << "\t";
      double df = fabs(funt - func);
      for (auto &i : direct)
      {
         stream << i << "\t";
      }
      stream << a << "\t";
      

      stream << fabs(cur[0] - tmp[0]) << "\t" << fabs(cur[1] - tmp[1]) << "\t";
      stream << df << "\t";
      stream << acos(scalar(direct,tmp)/(sqrt(scalar(direct,direct)) * sqrt(scalar(tmp, tmp)))) << "\t";
      
      for(auto &i : grad)
      {
         stream << i << "\t";
      }


      if (mod(cur,tmp) < eps)
      {
         break;
      }
      cur = tmp;
      grad = gf(cur);
      double ngg = scalar(grad,grad);
      double beta =  ngg/ gg;

      for(size_t i = 0; i < direct.size();i++)
      {
         direct[i] = -grad[i] + beta * direct[i];
      }
      if (fabs(scalar(grad, direct)) < 1e-12)
      {
         for (size_t i = 0; i < direct.size(); i++)
            direct[i] = -grad[i];
      }
      gg = ngg;
      k++;

      funt = func;
      stream << std::endl;
   }
   while (gg > eps* eps && maxI > k);

   for (auto &i : cur)
   {
      stream << std::setprecision(16) << i << "\t";
   }

   stream << std::setprecision(16) << funt << "\t";

   for (size_t i = 0; i < 7;i++)
   {
      stream << "-" << "\t";
   }

   for (auto &i : grad)
   {
      stream << std::setprecision(16) << i << "\t";
   }
   stream << std::endl;
   return cur;


}

void Ax(const std::vector<double> &A,const std::vector<double> &x,std::vector<double> &b)
{
   size_t size = x.size();
   b.resize(size);
   for(size_t i = 0; i < size; i++)
   {
      double buff = 0;
      for(size_t j = 0; j < size;j++)
      {
         buff += A[j + i * size] * x[j];
      }
      b[i] = buff;
   }



}

std::vector<double> Broiden
(
   std::function<double(const std::vector<double> &)> fun,
   std::function<std::vector<double>(const std::vector<double> &)> gf,
   const std::vector<double> &x0,
   double eps,double leps, size_t maxI,std::ostream &stream
)
{
   size_t size = x0.size();
   std::vector<double> B(size * size, 0);
   for(size_t i = 0; i < size;i++)
   {
      B[i + i * size] = 1;
   }
   std::vector<double> grad(size);
   grad = gf(x0);
   std::vector<double> direct(size, 0);
   std::vector<double> cur(x0);
   size_t k = 0;

   std::function<double(double)> layer = [&fun, &cur, &direct](double lambda)
      {
         std::vector<double> arg(cur);
         for (size_t i = 0; i < arg.size(); i++)
         {
            arg[i] += lambda * direct[i];
         }
         return fun(arg);
      };
   //stream << "x" << "\t" << "y" << "\t" << "f" << "\t" << "dirX" << "\t" << "dirY" << "\t" << "lambda" << "\t" << "dx" << "\t" << "dy" << "\t" << "df" << "\t" << "phi" << "\t" << "gradX" << "\t" << "gradY" << std::endl;
   while(scalar(grad,grad) > eps* eps && maxI > k)
   {
      for(auto &a : cur)
      {
         stream << std::setprecision(16) << a << "\t";
      }
      double func = fun(cur);
      stream << std::setprecision(16) << func << "\t";

      
      k++;

      Ax(B, grad, direct);
      for(auto &a : direct)
      {
         a = -a;
      }

      for (auto &a : direct)
      {
         stream << std::setprecision(16) << a << "\t";
      }

      std::pair<double, double> interval = findInterval(layer, 0.5);
      double a = squareInt(interval, leps, layer).second;
      stream << std::setprecision(16) << a << "\t";
      std::vector < double> s(size, 0);
      
      double phi = acos(scalar(direct, cur) / (sqrt(scalar(direct, direct)) * sqrt(scalar(cur, cur))));
      for(size_t i = 0; i < size;i++)
      {
         s[i] = direct[i] * a;
         cur[i] += s[i];
      }
      for (auto &i : s)
      {
         stream << std::setprecision(16) << i << "\t";
      }
      stream << std::setprecision(16) << fabs(func - fun(cur)) << "\t";

      stream << std::setprecision(16) << phi << "\t";
      for (auto &a : grad)
      {
         stream << std::setprecision(16) << a << "\t";
      }

      for (auto &a : B)
      {
         stream << std::setprecision(16) << a << "\t";
      }

      std::vector<double> y(grad);

      grad = gf(cur);
      for (size_t i = 0; i < size; i++)
      {
         y[i] = grad[i] - y[i];
      }
      std::vector<double> By(size, 0);
      Ax(B, y, By);
      
      std::vector<double> diff(size,0);
      for(size_t i = 0; i < size; i++)
      {
         diff[i] = s[i] - By[i];
      }
      double den = scalar(diff, y);


      for(size_t i = 0; i < size; i++)
      {
         for(size_t j = 0; j < size; j++)
         {
            B[i + j * size] += diff[i] * diff[j] / den;
         }
      }

      stream << std::endl;



   }

   for (auto &a : cur)
   {
      stream << std::setprecision(16) << a << "\t";
   }
   stream << std::setprecision(16) << fun(cur) << "\t";

   for (size_t i = 0; i < 7; i++)
   {
      stream << "-" << "\t";
   }

   for (auto &a : grad)
   {
      stream << std::setprecision(16) << a << "\t";
   }
   stream << std::endl;
   return cur;
}