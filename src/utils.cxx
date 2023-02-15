#include "utils.hpp"


std::string NWhiteSpaces(int n) { return std::string(std::max(0,n), ' ' );};
TString CompleteWhiteSpaces(TString text, int ntot) { return text+NWhiteSpaces(ntot-text.Length());};

bool FindInString(const std::string& search, const std::string& str) {return str.find(search)!=std::string::npos ;}

void PrintLoading(TString type, TString name, TString objname, std::string color) {
  Print("Loading "+type+" for "+CompleteWhiteSpaces(name)+": "+objname+"..."+reset, color);
}

void PrintParameter(TString name, double val, double err, std::string color){
  PrintLine("  --> "+CompleteWhiteSpaces(name)+Form(": %+5.2f +/- %5.2f", val, err), color);
}

void RemoveZerosFromGraph(TGraphErrors *graph) {
  for (int i = graph->GetN()-1; i != -1; --i) {// Must go backwards to keep point ordering
    if (graph->GetY()[i]==0 && graph->GetEY()[i]==0) {
      graph->RemovePoint(i);
    }
  }
}

void TruncateGraph(TGraphErrors *graph, const double& min, const double& max){
  for (int i = graph->GetN()-1; i != -1; --i) {// Must go backwards to keep point ordering
    double x = graph->GetX()[i];
    if (min>0 && x<min) graph->RemovePoint(i);
    if (max>0 && x>max) graph->RemovePoint(i);
  }
}


double oplus(double a, double b) {return sqrt(a*a + b*b);};


double fitError(Double_t *x, Double_t *p) {return -1000;}

std::function<double(Double_t *, Double_t *)> fitError_wrapper(TF1* func, TMatrixD err_matrix){
  return [func, err_matrix](Double_t *x, Double_t *p) -> double {

    int n_pars = func->GetNpar();
    // Partial derivatives as differentials with 10% step size
    std::vector<double> grad(n_pars);
    func->GradientPar(x, &grad[0]);
    // Perform standard error propagation
    double sumerr2(0);
    for (int i = 0; i != n_pars; ++i) {
      for (int j = 0; j != n_pars; ++j) {
        sumerr2 += err_matrix[i][j]*grad[i]*grad[j];
      }
    }
    return (func->Eval(*x) + p[0]*sqrt(sumerr2));
  };
}


void FuncToGraph(TF1* func, TMatrixD err_matrix, TGraphErrors *graph, double k){
  assert(func); assert(graph);
  std::unique_ptr<TGraph>gr; gr.reset(new TGraph(func));
  for (int i = 0; i < gr->GetN(); ++i) {
    double x = gr->GetX()[i];
    graph->SetPoint(i, x, gr->GetY()[i]);
    double err = fitError_wrapper(func, err_matrix)(&x,&k);
    graph->SetPointError(i, 0., err - gr->GetY()[i]);
  }
}

void FuncToHist(TF1* func, TMatrixD err_matrix, TH1D* hist, double k){
  assert(func); assert(hist);
  for (int j = 1; j != hist->GetNbinsX()+1; ++j) {
    double x = hist->GetBinCenter(j);
    hist->SetBinContent(j, x, func->Eval(x));
    double err = fitError_wrapper(func, err_matrix)(&x,&k);
    hist->SetBinError(j, 0., err - func->Eval(x));
  }
}


void multiplyGraph(TGraphErrors *graph, double scale) {
  for (int i = 0; i != graph->GetN(); ++i) {
    graph->SetPoint(i, graph->GetX()[i], scale*graph->GetY()[i]);
    graph->SetPointError(i, graph->GetEX()[i], scale*graph->GetEY()[i]);
  }
}


void multiplyGraph(TGraphErrors *graph, TF1 *func) {
  for (int i = 0; i != graph->GetN(); ++i) {
    double val = func->Eval(graph->GetX()[i]);
    graph->SetPoint(i, graph->GetX()[i], val*graph->GetY()[i]);
    graph->SetPointError(i, graph->GetEX()[i], val*graph->GetEY()[i]);
  }
}

void PropagateErrorToGraph(TGraphErrors *graph, std::map<int, TF1*> funcs, TMatrixD err_matrix) {
  for (int i = 0; i != graph->GetN(); ++i) {
    double x = graph->GetX()[i];
    double xerr = graph->GetEX()[i];

    double sumerr2(0);
    for (auto [i,func1]: funcs){
      for (auto [j,func2]: funcs){
        sumerr2 += func1->Eval(x)*func2->Eval(x)*err_matrix[i][j];
      }
    }
    graph->SetPointError(i, xerr, sqrt(sumerr2));
  }
}


TGraphErrors* MergeGraphs(TGraphErrors *graph1, TGraphErrors *graph2){
  std::vector<double> x_vals, y_vals, x_errs, y_errs;
  for (int bin = 0; bin < graph1->GetN(); bin++) {
    x_vals.push_back(graph1->GetX()[bin]);
    y_vals.push_back(graph1->GetY()[bin]);
    x_errs.push_back(graph1->GetEX()[bin]);
    y_errs.push_back(graph1->GetEY()[bin]);
  }
  if (!std::is_sorted(x_vals.begin(), x_vals.end())){
    throw std::invalid_argument("Unsorted graph used.");
  }
  if (x_vals.size()==0) return (TGraphErrors*)graph2->Clone();
  for (int bin = 0; bin < graph2->GetN(); bin++) {
    int index = closest(x_vals,graph2->GetX()[bin]);
    x_vals.insert(x_vals.begin()+index,graph2->GetX()[bin]);
    y_vals.insert(y_vals.begin()+index,graph2->GetY()[bin]);
    x_errs.insert(x_errs.begin()+index,graph2->GetEX()[bin]);
    y_errs.insert(y_errs.begin()+index,graph2->GetEY()[bin]);
  }
  TGraphErrors* out = new TGraphErrors(x_vals.size(), &x_vals[0], &y_vals[0], &x_errs[0], &y_errs[0]);
  return out;

}


bool IsParameterFixed (TF1* f1, int ipar) {
  double pmin,pmax;
  f1->GetParLimits(ipar,pmin,pmax);
  return (pmin*pmax !=0 && pmin >= pmax);
}
