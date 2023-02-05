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


PFCompositionContainer::PFCompositionContainer(TString name_, TString type_){
  set_name(name_);
  set_type(type_);
  raw_.reset(new TGraphErrors());
  input_.reset(new TGraphErrors());
  output_.reset(new TGraphErrors());
  variation_.reset(new TGraphErrors());
}


void PFCompositionContainer::add_raw(TGraphErrors* graph_, TString name_){
  raws_.emplace_back((TGraphErrors*)graph_->Clone(name_+"_raw"));
}

void PFCompositionContainer::add_graph(TGraphErrors* graph_, TString name_){
  inputs_.emplace_back((TGraphErrors*)graph_->Clone(name_+"_input"));
  outputs_.emplace_back((TGraphErrors*)graph_->Clone(name_+"_output"));
  variations_.emplace_back((TGraphErrors*)graph_->Clone(name_+"_variation"));
}

void PFCompositionContainer::set_combination(){
  for (auto &x: raws_)       raw_.reset(MergeGraphs(raw_.get(), x.get()));
  for (auto &x: inputs_)     input_.reset(MergeGraphs(input_.get(), x.get()));
  for (auto &x: outputs_)    output_.reset(MergeGraphs(output_.get(), x.get()));
  for (auto &x: variations_) variation_.reset(MergeGraphs(variation_.get(), x.get()));
  raw_->SetName(type()+"_raw");
  input_->SetName(type()+"_input");
  output_->SetName(type()+"_output");
  variation_->SetName(type()+"_variation");

}

ostream& operator<<(ostream& os, const PFCompositionContainer& pf) {
  PrintLine("name: "+pf.name(), pf.color());
  PrintLine("  --> type: "+pf.type(), pf.color());
  Print("  --> # inputs: "+to_string(pf.raws().size())+reset, pf.color());
  return os;
}

ShapeContainer::ShapeContainer(TString name_, TString form_, TString appliesTo_, int index_, bool ispositive_){
  set_name(name_);
  set_form(form_);
  set_appliesTo(appliesTo_);
  set_index(index_);
  set_ispositive(ispositive_);
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


SystematicContainer::SystematicContainer(TString name_, TString appliesTo_, int index_, TString hname_, TH1D* hist_){
  set_name(name_);
  set_appliesTo(appliesTo_);
  set_index(index_);
  set_hname(hname_);
  set_hist((TH1D*)(hist_)->Clone(name_));
}

ostream& operator<<(ostream& os, const SystematicContainer& sys) {
  PrintLine("name: "+sys.name(), sys.color());
  PrintLine("  --> appliesTo: "+sys.appliesTo(), sys.color());
  PrintLine("  --> hname: "+sys.hname(), sys.color());
  Print("  --> index: "+to_string(sys.index())+reset, sys.color());
  return os;
}
