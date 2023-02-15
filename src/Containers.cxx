#include "Containers.hpp"
#include "constants.hpp"

using namespace std;

DataContainer::DataContainer(TString name_, TString type_, TString hname_, TGraphErrors* graph_){
  set_name(name_);
  set_type(type_);
  set_hname(hname_);
  RemoveZerosFromGraph(graph_);
  // TruncateGraph(graph_,15);
  set_raw((TGraphErrors*)graph_->Clone(name_+"_raw"));
  set_input((TGraphErrors*)graph_->Clone(name_+"_in"));
  set_output((TGraphErrors*)graph_->Clone(name_+"_out"));
  set_variation((TGraphErrors*)graph_->Clone(name_+"_variation"));
}

ostream& operator<<(ostream& os, const DataContainer& dt) {
  PrintLine("name: "+dt.name(), dt.color());
  PrintLine("  --> type: "+dt.type(), dt.color());
  Print("  --> hname: "+dt.hname()+reset, dt.color());
  return os;
}


ShapeContainer::ShapeContainer(TString name_, TString form_, TString appliesTo_, int index_, bool ispositive_, bool freeze_, double initial_){
  set_name(name_);
  set_form(form_);
  set_appliesTo(appliesTo_);
  set_index(index_);
  set_freeze(freeze_);
  set_initial(initial_);
  set_func(new TF1("f1_"+name_+"_"+appliesTo_, form_, func_range_min,func_range_max));
}

ostream& operator<<(ostream& os, const ShapeContainer& fitshape) {
  PrintLine("name: "+fitshape.name(), fitshape.color());
  PrintLine("  --> form: "+fitshape.form(), fitshape.color());
  PrintLine("  --> appliesTo: "+fitshape.appliesTo(), fitshape.color());
  PrintLine("  --> index: "+to_string(fitshape.index()), fitshape.color());
  Print("  --> ispositive: "+to_string(fitshape.ispositive())+reset, fitshape.color());
  return os;
}


NuisanceContainer::NuisanceContainer(TString name_, TString appliesTo_, int index_, TString hname_, TH1D* hist_){
  set_name(name_);
  set_appliesTo(appliesTo_);
  set_index(index_);
  set_hname(hname_);
  set_hist((TH1D*)(hist_)->Clone(name_));
}

ostream& operator<<(ostream& os, const NuisanceContainer& sys) {
  PrintLine("name: "+sys.name(), sys.color());
  PrintLine("  --> appliesTo: "+sys.appliesTo(), sys.color());
  PrintLine("  --> hname: "+sys.hname(), sys.color());
  Print("  --> index: "+to_string(sys.index())+reset, sys.color());
  return os;
}
