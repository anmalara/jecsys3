#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>

#include "TF1.h"
#include "TH1D.h"
#include "TMatrixD.h"
#include "TGraphErrors.h"

const std::string red("\x1b[0;31m");
const std::string green("\x1b[0;32m");
const std::string yellow("\x1b[0;33m");
const std::string blue("\x1b[0;34m");
const std::string magenta("\x1b[0;35m");
const std::string cyan("\x1b[0;36m");
const std::string bold("\033[1m");
const std::string reset("\x1b[0m");


std::string NWhiteSpaces(int n);
TString CompleteWhiteSpaces(TString text, int ntot=20);

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
  for (const auto& x: vec) { os << x << " ";}
  return os;
}

template<class T>
void Print(const T text, const std::string color){
  std::cout << color << text << " ";
}

template<class T>
void Print(const T text){
  std::cout << text << " ";
}

template<class T>
void PrintLine(const T text, const std::string color){
  Print(text,color); std::cout << reset << std::endl;
}

template<class T>
void PrintLine(const T text){
  Print(text); std::cout << reset << std::endl;
}

void PrintLoading(TString type, TString name, TString objname, std::string color=yellow);
void PrintParameter(TString name, double val, double err, std::string color=blue);

bool FindInString(const std::string& search, const std::string& str);

template <typename T>
int FindInVector(const std::vector<T>& vec, const T& el) {
  int index = -1;
  // Find given element in vector
  auto it = std::find(vec.begin(), vec.end(), el);
  if (it != vec.end()) index = distance(vec.begin(), it);
  return index;
}

template <typename T>
bool FindInMap(const std::map<TString,T>& map, const TString& el) {
  return (map.find(el) != map.end());
}

#include <sstream>

template <typename T>
std::string to_string_with_precision(const T value, const int n = 3) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << value;
    return out.str();
}


void RemoveZerosFromGraph(TGraphErrors *graph);
void TruncateGraph(TGraphErrors *graph, const double& min = -1, const double& max=-1);

double oplus(double a, double b);

double fitError(Double_t *x, Double_t *p);
std::function<double(Double_t *, Double_t *)> fitError_wrapper(TF1* func, TMatrixD err_matrix);

void FuncToGraph(TF1* func, TMatrixD err_matrix, TGraphErrors* graph, double k_err=1);
void FuncToHist(TF1* func, TMatrixD err_matrix, TH1D* hist, double k_err=1);

void multiplyGraph(TGraphErrors *graph, double scale);
void multiplyGraph(TGraphErrors *graph, TF1 *func);
void PropagateErrorToGraph(TGraphErrors *graph, std::map<int, TF1*> funcs, TMatrixD err_matrix);
TGraphErrors* MergeGraphs(TGraphErrors *graph1, TGraphErrors *graph2);

template <typename T>
int closest(std::vector<T> const& vec, T value) {
    auto const it = std::lower_bound(vec.begin(), vec.end(), value);
    T val = (it == vec.end())? -1 : *it;
    return find(vec.begin(), vec.end(),val) - vec.begin();
}


bool IsParameterFixed (TF1* f1, int ipar);
