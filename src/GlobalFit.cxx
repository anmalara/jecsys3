#include "GlobalFit.hpp"
#include "Containers.hpp"
#include "utils.hpp"

using namespace std;

// Set to a dummy initial value
TString GlobalFit::current_obs = "";

GlobalFit::GlobalFit(string input_json) {
  // Constructor: read info from json

  ifstream f_(input_json);
  nlohmann::json infos = nlohmann::json::parse(f_);
  runName = infos["run"];
  mode = infos["mode"];
  eta_min = infos["eta_min"];
  eta_max = infos["eta_max"];
  if (infos.contains("output_fname"))
    output_fname = infos["output_fname"];

  for (auto i : infos["samples"].get<vector<string>>())
    samples.push_back(i);
  for (auto i : infos["hdm_methods"].get<vector<string>>())
    hdm_methods.push_back(i);
  for (auto i : infos["types"].get<vector<string>>())
    types.push_back(i);
  for (auto i : infos["shapes"].get<vector<string>>())
    shapes_allowed.push_back(i);

  f_.close();

  PrintLine("Running on:", green);
  PrintLine("  --> run :" + runName, green);
  PrintLine("  --> mode :" + mode, green);
  PrintLine("  --> eta_min :" + eta_min, green);
  PrintLine("  --> eta_max :" + eta_max, green);
  Print("  --> samples:", green);
  PrintLine(samples);
  Print("  --> hdm_methods:", green);
  PrintLine(hdm_methods);
  Print("  --> types:", green);
  PrintLine(types);
  Print("  --> shapes:", green);
  PrintLine(shapes_allowed);
  for (auto type : types) {
    input_hnames.push_back(type);
    for (auto sample : samples) {
      for (auto method : hdm_methods) {
        input_hnames.push_back(type + "_" + sample + "_" + method);
      }
    }
  }

  OpenFiles();
}

GlobalFit::~GlobalFit() {
  // Destructor: takes care of closing files
  for (auto [fname, f_] : input_files) {
    if (debug)
      PrintLine("Closing file for " + fname, yellow);
    f_->Close();
  }
  for (auto [fname, f_] : output_files) {
    if (debug)
      PrintLine("Closing file for " + fname, yellow);
    f_->Close();
  }
}

void GlobalFit::OpenFiles() {
  // Open input/output files
  for (auto [name, fname] : input_fnames) {
    if (debug)
      PrintLoading("file", name, ReplaceDefault(fname));
    input_files[name] = new TFile(ReplaceDefault(fname), "READ");
    assert(input_files[name] && !input_files[name]->IsZombie());
    if (debug)
      PrintLine("successfully", yellow);
  }
  output_files["output"] = new TFile(ReplaceDefault(output_fname), "RECREATE");
}

TString GlobalFit::ReplaceDefault(TString name) {
  // Function to update the default strings
  TString res = name;
  res = res.ReplaceAll("RUN", runName).ReplaceAll("MODE", mode);
  res = res.ReplaceAll("ETAMIN", eta_min).ReplaceAll("ETAMAX", eta_max);
  return res;
};

void GlobalFit::LoadInputs() {
  // Load data points

  for (auto [name, info] : input_hnames_map) {
    if (FindInVector(input_hnames, name) < 0) {
      if (debug)
        PrintLine("Skipping hist: " + name, yellow);
      continue;
    }
    TString hname = ReplaceDefault(info["hname"]);
    if (debug)
      PrintLoading("hist", name, hname);
    unique_ptr<TGraphErrors> g;
    g.reset((TGraphErrors *)input_files[info["fname"].Data()]->Get(hname));
    my_data[name] = new DataContainer(name, info["type"], hname, g.get());
    nTotalPoints += my_data[name]->GetN();
    if (debug) {
      PrintLine("successfully", yellow);
      cout << *my_data[name] << endl;
    }
  }

  // Special treatment of recoild: TODO
  for (auto [name, info] : recoil_hnames_map) {
    TString hname = ReplaceDefault(info["hname"]);
    if (debug)
      PrintLoading("recoil", name, hname);
    unique_ptr<TGraphErrors> g;
    g.reset((TGraphErrors *)input_files[info["fname"].Data()]->Get(hname));
    recoils[name] = new DataContainer(name, info["type"], hname, g.get());
    if (debug) {
      PrintLine("successfully", yellow);
      cout << *recoils[name] << endl;
    }
  }
}

void GlobalFit::LoadSystematics() {
  // Load Nuisances
  for (auto [name, info] : nuisances_map) {
    TString appliesTo = info["appliesTo"];
    if (FindInVector(input_hnames, appliesTo) < 0) {
      if (debug)
        PrintLine("Skipping nuisance: " + name + " index: " + to_string(FindInVector(input_hnames, appliesTo)), yellow);
      continue;
    }
    TString hname = ReplaceDefault(info["hname"]);
    if (debug)
      PrintLoading("nuisance", name, hname);
    unique_ptr<TH1D> h;

    if (FindInString("scale", name.Data())) {
      h.reset(new TH1D(name, name, 10000, 10, 5000));
      double scale = stod(hname.ReplaceAll("scale", "").Data());
      hname = name;
      for (int i = 1; i != h->GetNbinsX() + 1; ++i) {
        h->SetBinContent(i, scale);
        h->SetBinError(i, scale);
      }
    } else {
      h.reset((TH1D *)input_files[info["fname"].Data()]->Get(hname));
    }
    nuisances[name] = new NuisanceContainer(name, appliesTo, nuisances.size(), hname, h.get());
    assert(nuisances[name]);
    if (debug) {
      PrintLine("successfully", yellow);
      cout << *nuisances[name] << endl;
    };
  }
}

void GlobalFit::LoadShapes() {
  // Load functional forms for POIs
  for (auto [name, info] : shapes_map) {
    TString appliesTo = info["appliesTo"];
    if (FindInVector(types, appliesTo) < 0 || FindInVector(shapes_allowed, TString(name).ReplaceAll("_" + appliesTo, "")) < 0) {
      if (debug)
        PrintLine("Skipping shape: " + name + " index: " + to_string(FindInVector(types, appliesTo)) + " " + to_string(FindInVector(shapes_allowed, appliesTo)),
                  yellow);
      continue;
    }
    if (debug)
      PrintLoading("shape", name, "");
    int index = (info["appliesTo"] == "Resp") ? FindInVector(shapes_allowed, name) : shapes[info["type"]]->index();
    shapes[name] =
        new ShapeContainer(name, info["form"], info["appliesTo"], index, atoi(info["ispositive"]), atoi(info["freeze"]), stod(info["initial"].Data()));
    assert(shapes[name]);
    if (debug) {
      PrintLine("successfully", yellow);
      cout << *shapes[name] << endl;
    };
    shape_types.insert(info["type"]);
  }
}

void GlobalFit::LoadReference() {
  for (auto [name, info] : reference_obj_map) {
    TString hname = ReplaceDefault(info["hname"]);
    if (debug)
      PrintLoading("reference", name, hname);
    reference_objects[name] = (TH1D *)input_files[info["fname"].Data()]->Get(hname)->Clone(name);
    assert(reference_objects[name]);
    reference_objects[name]->SetDirectory(0);
    if (debug) {
      PrintLine("successfully", yellow);
    }
  }
}

void GlobalFit::LoadFSR() {
  // Load fsr corrections
  for (auto [name, info] : kfsr_hnames_map) {
    TString appliesTo = info["appliesTo"];
    if (FindInVector(input_hnames, appliesTo) < 0) {
      if (debug)
        PrintLine("Skipping fsr: " + name + " index: " + to_string(FindInVector(input_hnames, appliesTo)), yellow);
      continue;
    }
    TString hname = ReplaceDefault(info["hname"]);
    if (debug)
      PrintLoading("kfsr", name, hname);
    unique_ptr<TH1D> h;
    h.reset((TH1D *)input_files[info["fname"].Data()]->Get(hname));
    fsrs[name] = new NuisanceContainer(name, appliesTo, -1, hname, h.get());
    assert(fsrs[name]);
    if (debug) {
      PrintLine("successfully", yellow);
      cout << *fsrs[name] << endl;
    };
  }
}

void GlobalFit::ScaleFSR() {
  // Scale raw with frs corrections (HDM method)
  for (auto [name, fsr] : fsrs) {
    auto dt = my_data[fsr->appliesTo()];
    TGraphErrors *raw = dt->raw();
    TGraphErrors *input = dt->input();
    TGraphErrors *output = dt->output();
    TGraphErrors *variation = dt->variation();
    for (int i = 0; i != raw->GetN(); ++i) {
      double pt = raw->GetX()[i];
      double r = raw->GetY()[i];
      double shift = fsr->hist()->GetBinContent(fsr->hist()->FindBin(pt));

      if (FindInString("multijet", name.Data())) {
        TString recoil_name = dt->name().ReplaceAll("multijet", "multijet_recoil");
        TGraphErrors *recoil_raw = recoils[recoil_name]->raw();
        TGraphErrors *recoil_input = recoils[recoil_name]->input();
        TGraphErrors *recoil_output = recoils[recoil_name]->output();
        TGraphErrors *recoil_variation = recoils[recoil_name]->variation();
        TH1D *herr = reference_objects["herr"];
        double err_r = raw->GetEX()[i];
        double err_pt = raw->GetEY()[i];

        double ptref = recoil_raw->GetY()[i] * pt;
        double jesref = herr->GetBinContent(herr->FindBin(ptref));
        double jes = herr->GetBinContent(herr->FindBin(pt));
        cout << "MJS " << pt << " " << ptref << " " << recoil_raw->GetX()[i] << " " << jesref << " " << jes << " ";
        cout << r << " " << shift << " " << jesref * (r + shift) << " " << (jesref * r) + shift << endl;
        raw->SetPoint(i, pt, jesref * r);
        input->SetPoint(i, pt, jesref * (r + shift));
        output->SetPoint(i, pt, jesref * (r + shift));
        variation->SetPoint(i, pt, jesref * (r + shift));

        recoil_raw->SetPoint(i, ptref, jes / r);
        recoil_input->SetPoint(i, ptref, jes / (r + shift));
        recoil_output->SetPoint(i, ptref, jes / (r + shift));
        recoil_variation->SetPoint(i, ptref, jes / (r + shift));
        recoil_raw->SetPointError(i, err_r, err_pt);
        recoil_input->SetPointError(i, err_r, err_pt);
        recoil_output->SetPointError(i, err_r, err_pt);
        recoil_variation->SetPointError(i, err_r, err_pt);
      } else {
        input->SetPoint(i, pt, r + shift);
        output->SetPoint(i, pt, r + shift);
        variation->SetPoint(i, pt, r + shift);
      }
    }
  }
}

void GlobalFit::CleanGraphs() {
  // Helper function to remove points from graphs
  for (auto [name, recoil] : recoils) {
    TruncateGraph(recoil->input(), ptmin_multijet, ptmax_multijet);
    TruncateGraph(recoil->output(), ptmin_multijet, ptmax_multijet);
    TruncateGraph(recoil->variation(), ptmin_multijet, ptmax_multijet);
  }

  for (auto [dt_name, dt] : my_data) {
    if (dt->type() == "Resp")
      continue;
    int removed_points = 0;
    for (int i = 0; i < dt->input()->GetN(); ++i) {
      if (dt->input()->GetX()[i] <= ptmin_pf)
        removed_points++;
      if (dt->input()->GetX()[i] >= ptmax_pf)
        removed_points++;
    }
    nTotalPoints -= removed_points;
    TruncateGraph(dt->input(), ptmin_pf, ptmax_pf);
    TruncateGraph(dt->output(), ptmin_pf, ptmax_pf);
    TruncateGraph(dt->variation(), ptmin_pf, ptmax_pf);
  }
}

void GlobalFit::SetupFitFunction() {
  // Set the minimizer tool for the global fit
  if (debug)
    PrintLine("Loaded " + to_string(my_data.size()) + " graphs", green);
  if (debug)
    PrintLine("Loaded " + to_string(shapes.size()) + " shapes", green);
  if (debug)
    PrintLine("Loaded " + to_string(nuisances.size()) + " nuisances", green);

  nFitPars = shape_types.size();
  nNuisancePars = nuisances.size();
  nTotPars = nFitPars + nNuisancePars;

  auto jesFit_wrapper_ = jesFit_wrapper(reference_objects["hjesref"], shapes);
  _jesFit = new TF1("jesFit", jesFit_wrapper_, fit_min, fit_max, nFitPars);
  for (auto [name, shape] : shapes) {
    int index = shape->index();
    if (index > nFitPars)
      continue;
    if (shape->freeze()) {
      _jesFit->FixParameter(index, shape->initial());
      nFrozenPars++;
    } else
      _jesFit->SetParameter(index, shape->initial());
  }

  // TMinuit *fitter = new TMinuit(nTotPars);
  fitter = new TFitter(nTotPars);
  fitter->SetFCN(jesFitter);
  for (int i = 0; i != nTotPars; ++i) {
    fitter->SetParameter(i, "", _jesFit->GetParameter(i), (i < nFitPars ? ScaleFullSimShape : 1), -100, 100);
    if (i < nFitPars && IsParameterFixed(_jesFit, i))
      fitter->FixParameter(i);
  }

  PrintLine("Global fit has " + to_string(nTotPars) + " total parameters:", blue);
  PrintLine("  --> " + to_string(nFitPars) + " fit parameters", blue);
  PrintLine("    --> " + to_string(nFrozenPars) + " parameters are frozen", blue);
  PrintLine("  --> " + to_string(nNuisancePars) + " nuisance parameters", blue);
  if (penalizeFitPars)
    PrintLine("  --> Fit parameters have Gaussian prior", blue);
  PrintLine("  --> " + to_string(nTotalPoints) + " data points", blue);
}

void GlobalFit::DoGlobalFit() {

  // Run fitter (multiple times if needed)
  for (int i = 0; i != number_of_fit_iterations; ++i) {
    fitter->ExecuteCommand("MINI", 0, 0);
  }

  // Verify that the degrees of freedom make sense. Important to check
  // immediately.
  nFittedDataPoints -= nNuisancePars;
  if (penalizeFitPars)
    nFittedDataPoints -= nFitPars - nFrozenPars;
  assert(nFittedDataPoints == nTotalPoints);

  // Set the error matrix
  error_matrix.reset(new TMatrixD(nTotPars, nTotPars));
  gMinuit->mnemat(error_matrix->GetMatrixArray(), nTotPars);

  // Retrieve the chi2 for the individual components
  Double_t tmp_par[nTotPars], grad[nTotPars];
  Double_t chi2_gbl(0);
  for (int i = 0; i < nTotPars; ++i)
    tmp_par[i] = fitter->GetParameter(i);
  jesFitter(nTotPars, grad, chi2_gbl, tmp_par, 1);

  // Retrieve the chi2 for the individual components
  Double_t chi2_nuis(0), chi2_par(0), chi2_data(0);
  int npar_true(0), nnuis_true(0);

  for (int i = 0; i != nTotPars; ++i) {
    double val = fitter->GetParameter(i);
    double err = fitter->GetParError(i);
    if (fabs(val) != 0 || fabs(err - 1) > 1e-2) {
      if (i < nFitPars) {
        ++npar_true;
        chi2_par += pow(val, 2);
      } else {
        ++nnuis_true;
        chi2_nuis += pow(val, 2);
      }
    }
  }

  // Verify that numbers of parameters make sense
  assert(nFitPars == npar_true);
  assert(nNuisancePars == nnuis_true);

  for (auto [dt_name, dt] : my_data) {
    TGraphErrors *graph_output = dt->output();
    assert(graph_output);
    current_obs = dt->type(); // this is needed inside _jesFit
    for (int j = 0; j != graph_output->GetN(); ++j) {
      double x = graph_output->GetX()[j];
      double y = graph_output->GetY()[j];
      double ey = graph_output->GetEY()[j];
      chi2_data += pow((y - _jesFit->Eval(x)) / ey, 2);
    }
  }

  // Summary of the global fit
  PrintLine("Output GlobalFit", blue);
  PrintLine("  --> " + to_string(nTotalPoints) + " data points used", blue);
  PrintLine("  --> " + to_string(nFitPars) + " fit parameters", blue);
  PrintLine("    --> " + to_string(nFrozenPars) + " parameters are frozen", blue);
  PrintLine("  --> " + to_string(nnuis_true) + " nuisances", blue);
  PrintLine(Form("  --> fitting range [%1.0f,%1.0f]", _jesFit->GetXmin(), _jesFit->GetXmax()), blue);
  PrintLine(Form("  --> Total     chi2/NDF  = %1.1f / %d", chi2_gbl, nFittedDataPoints), blue);
  PrintLine(Form("  --> Data      chi2/NDF  = %1.1f / %d", chi2_data, nTotalPoints), blue);
  PrintLine(Form("  --> Nuisance  chi2/Nnuis = %1.1f / %d", chi2_nuis, nNuisancePars), blue);
  PrintLine(Form("  --> Parameter chi2/Npar = %1.1f / %d", chi2_par, nFitPars), blue);

  // Verify that chi2 and d.o.f. make sense
  //   assert(std::fabs(chi2_gbl - (chi2_data + chi2_nuis + chi2_par)) < 1e-6); TODO. this must make sense also with parameters
  assert(nFittedDataPoints == nTotalPoints + nFitPars + nNuisancePars);

  PrintLine("Fitted parameters (for Resp):", blue);
  for (auto [name, shape] : shapes) {
    if (shape->appliesTo() != "Resp")
      continue;
    double val = fitter->GetParameter(shape->index());
    double err = fitter->GetParError(shape->index());
    PrintParameter(name, val, err);
  }

  PrintLine("Nuisance parameters:", blue);
  for (auto [name, nuis] : nuisances) {
    double val = fitter->GetParameter(nFitPars + nuis->index());
    double err = fitter->GetParError(nFitPars + nuis->index());
    PrintParameter(name, val, err);
  }

  PrintLine("Error matrix:", blue);

  for (int i = 0; i != nFitPars; ++i) {
    for (int j = 0; j != nFitPars; ++j) {
      double err = (*error_matrix.get())[i][j];
      Print(Form((err >= 0) ? " +%2.4f" : " %2.4f", err), blue);
    }
    PrintLine("");
  }
  PrintLine("");
}

void GlobalFit::StoreFitOutput() {
  // Store input and output into file
  output_files["output"]->cd();

  // Set range and granularity for plotting
  _jesFit->SetRange(func_range_min, func_range_max);
  _jesFit->SetNpx(func_range_max - func_range_min);
  auto fitError_wrapper_ = fitError_wrapper(_jesFit, *error_matrix.get());
  TF1 *fit_err = new TF1("fitError", fitError_wrapper_, func_range_min, func_range_max, 1);
  fit_err->SetNpx(func_range_max - func_range_min);

  // Convert global fit function into hist and graph
  // clang-format off
  double binning[] = {/*1, 5, 6, 8,*/ 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133,
  153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
  1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941,
  3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
  // clang-format on
  const double nbins = sizeof(binning) / sizeof(binning[0]) - 1;

  for (auto type : types) {
    current_obs = type; // this is needed inside _jesFit
    // double scale = (type == "Resp") ? 1. : 1. / ScaleFullSimShape;
    double scale = 1;
    // convert TF1 to TGraphErrors
    TGraphErrors *graph = new TGraphErrors();
    FuncToGraph(_jesFit, *error_matrix.get(), graph);
    if (type == "Resp") {
      // TODO YOU'RE HERE. not sure why this is chaning the stored out func.
      for (int i = graph->GetN() - 1; i != -1; --i) {
        graph->SetPoint(i, graph->GetX()[i], graph->GetY()[i] - 1);
      }
      //   multiplyGraph(graph, 100);
    }
    multiplyGraph(graph, scale);
    // convert TF1 to TGraphErrors
    TH1D *hist = new TH1D("jesFit_hist_" + type, ";p_{T} (GeV);" + type, nbins, binning);
    hist->SetDirectory(0);
    FuncToHist(_jesFit, *error_matrix.get(), hist);
    hist->Scale(scale);
    // Store
    graph->Write("jesFit_graph_" + type, TObject::kOverwrite);
    hist->Write(hist->GetName(), TObject::kOverwrite);
    _jesFit->Write("jesFit_" + type, TObject::kOverwrite);
    fit_err->SetParameter(0, -1);
    fit_err->Write("jesFit_down_" + type, TObject::kOverwrite);
    fit_err->SetParameter(0, +1);
    fit_err->Write("jesFit_up_" + type, TObject::kOverwrite);
  }
  error_matrix->Write("error_matrix", TObject::kOverwrite);

  // Store all jes response
  std::vector<double> bins;
  for (auto [name, dt] : my_data) {
    TGraphErrors *raw = dt->raw();
    TGraphErrors *input = dt->input();
    TGraphErrors *output = dt->output();
    TGraphErrors *variation = dt->variation();

    for (int i = 0; i < output->GetN(); ++i)
      bins.push_back((int)output->GetX()[i]);
    double scale = (dt->type() == "Resp") ? 1. : 1. / ScaleFullSimShape; // transform PF comp in percentage
    multiplyGraph(input, scale);
    multiplyGraph(output, scale);
    multiplyGraph(variation, scale);

    std::map<int, TF1 *> funcs;
    for (auto [name_, shape] : shapes) {
      if (shape->appliesTo() == dt->type()) {
        funcs[shape->index()] = shape->func();
      }
    }
    raw->Write(name + "_raw", TObject::kOverwrite);
    input->Write(name + "_prefit", TObject::kOverwrite);
    output->Write(name + "_postfit", TObject::kOverwrite);
    variation->Write(name + "_variation_input", TObject::kOverwrite);
    PropagateErrorToGraph(variation, funcs, *error_matrix.get());
    variation->Write(name + "_variation_output", TObject::kOverwrite);
  }

  sort(bins.begin(), bins.end());
  bins.erase(unique(bins.begin(), bins.end()), bins.end());

  TGraphErrors *Resp_comb = new TGraphErrors();
  TGraphErrors *Resp_comb_shift = new TGraphErrors();
  for (int i = 0; i < bins.size(); i++) {
    double sum = 0, sum_shift = 0, sum_weight = 0;
    double x = bins[i];
    for (auto [name, dt] : my_data) {
      if (dt->type() != "Resp")
        continue;
      TGraphErrors *ref = dt->output();
      TGraphErrors *ref_shift = dt->input();
      if (x < ref->GetX()[0])
        continue;
      if (x > ref->GetX()[ref->GetN() - 1])
        continue;
      int i_graph = -1;
      for (int ii = 0; ii < ref->GetN() - 1; ii++) {
        if (x <= ref->GetX()[ii + 1]) {
          i_graph = ii;
          break;
        }
      }
      if (i_graph < 0)
        continue;
      double w = 1. / ref->GetEY()[i_graph];
      sum += w * ref_shift->Eval(x);
      sum_shift += w * ref->Eval(x);
      sum_weight += w;
    }
    if (sum_weight == 0 && sum == 0) {
      sum = sum_weight = 1;
    }
    Resp_comb->SetPoint(i, x, sum / sum_weight);
    Resp_comb_shift->SetPoint(i, x, (sum_shift / sum_weight - 1) / ScaleFullSimShape);

    Resp_comb->SetPointError(i, 0, 1 / sum_weight);
    Resp_comb_shift->SetPointError(i, 0, (1 / sum_weight) / ScaleFullSimShape);
  }

  Resp_comb->Write("Resp_comb", TObject::kOverwrite);
  Resp_comb_shift->Write("Resp_comb_shift", TObject::kOverwrite);

  for (auto [name, recoil] : recoils) {
    recoil->raw()->Write(name + "_raw", TObject::kOverwrite);
    recoil->input()->Write(name + "_prefit", TObject::kOverwrite);
    recoil->output()->Write(name + "_postfit", TObject::kOverwrite);
    recoil->variation()->Write(name + "_variation", TObject::kOverwrite);
  }

  // Store all nuisances
  for (auto [name, nuisance] : nuisances) {
    nuisance->hist()->SetYTitle(name);
    nuisance->hist()->Write("nuisance_" + name, TObject::kOverwrite);
  }

  // Store all shapes: pre/post fit
  for (auto [name, shape] : shapes) {
    current_obs = shape->appliesTo();
    shape->func()->Write("shape_input_" + name, TObject::kOverwrite);
    TF1 *prefit = new TF1("shape_prefit_" + name, "[0]*(" + shape->form() + ")", func_range_min, func_range_max);
    TF1 *postfit = new TF1("shape_postfit_" + name, "[0]*(" + shape->form() + ")", func_range_min, func_range_max);
    prefit->SetParameter(0, 1);
    postfit->SetParameter(0, _jesFit->GetParameter(shape->index()));
    prefit->Write(prefit->GetName(), TObject::kOverwrite);
    postfit->Write(postfit->GetName(), TObject::kOverwrite);
  }

  // Store reference objects
  for (auto [name, obj] : reference_objects) {
    obj->Write(name, TObject::kOverwrite);
  }
}

void GlobalFit::Run() {
  LoadInputs();
  LoadSystematics();
  LoadShapes();
  LoadReference();
  LoadFSR();
  ScaleFSR();
  CleanGraphs();
  SetupFitFunction();
  DoGlobalFit();
  StoreFitOutput();
}

function<double(Double_t *, Double_t *)> jesFit_wrapper(TH1D *hjesref, map<TString, ShapeContainer *> shapes) {

  return [hjesref, shapes](Double_t *x, Double_t *p) -> double {
    double var = (GlobalFit::current_obs == "Resp" ? 1. : 0.);
    double scale = (GlobalFit::current_obs == "Resp") ? 1. : 1. / GlobalFit::ScaleFullSimShape;
    double pt = x[0];

    double jesref = 1;
    if (GlobalFit::useJESref && GlobalFit::current_obs == "Resp") {
      assert(hjesref);
      jesref = hjesref->Interpolate(pt);
    }

    // // Load available shapes for this observable
    for (auto [name, shape] : shapes) {
      if (shape->appliesTo() != GlobalFit::current_obs)
        continue;
      // Calculate variation and add it to total
      int index = shape->index();
      assert(index >= 0);
      TF1 *f1 = shape->func();
      assert(f1);
      double par = p[index];
      if (shape->ispositive())
        par = max(par, 0.);
      // if (shape->ispositive()) par = min(1.,max(par,0.));
      var += par * f1->Eval(pt) * scale;
    return (var / jesref);
  };
}

// Dummy value: not used at the moment but it can be changes if needed
Double_t jesFit(Double_t *x, Double_t *p) { return -1000; }

void jesFitter(Int_t &npar, Double_t *grad, Double_t &chi2, Double_t *par, Int_t flag) {

  // Basic checks
  assert(_jesFit);

  // Parametes for nuisances (nuisances)
  int _nFitPars = _jesFit->GetNpar();
  int _nNuisancePars = npar - _nFitPars;
  Double_t *fit_pars = &par[0];
  Double_t *nuisance_pars = &par[_nFitPars];

  if (flag) {
    // do the following calculation:
    // chi2 = sum_i( (x_i+sum_s(a_s y_si) -fit)^2/sigma_i^2) + sum_s(a_s^2)
    chi2 = 0;
    nFittedDataPoints = 0;

    // Loop over input data (graphs x bins)
    // - for each point, add up source eigenvectors x nuisance parameters
    // - then calculate chi2 adding up residuals + nuisance parameters
    for (auto [dt_name, dt] : my_data) {
      // if (FindInString("recoil",dt_name.Data())) continue;
      // cout << *dt << endl;
      TString dt_type = dt->type();
      TGraphErrors *graph_input = dt->input();
      assert(graph_input);
      TGraphErrors *graph_output = dt->output();
      assert(graph_output);
      TGraphErrors *graph_variation = dt->variation();
      assert(graph_variation);

      DataContainer *dt_recoil;
      bool is_multijet_resp = dt_name.Contains("multijet") && dt_type == "Resp";
      if (is_multijet_resp) {
        dt_recoil = recoils[TString(dt_name).ReplaceAll("multijet", "multijet_recoil")];
      }

      for (int bin = 0; bin < graph_input->GetN(); ++bin) {

        // Retrieve central value and uncertainty for this point
        double pt = graph_input->GetX()[bin];
        double data_point = graph_input->GetY()[bin];
        double sigma = graph_input->GetEY()[bin];

        // Calculate fit value at this point
        GlobalFit::current_obs = dt_type;
        _jesFit->SetParameters(fit_pars);
        double fit = _jesFit->EvalPar(&pt, fit_pars);

        // // For multijet balancing, multiply data by reference JES
        // if (is_multijet_resp) {
        //   // what is happening?? Should be removed
        //   double ptref = dt_recoil->raw()->GetX()[bin];
        //   double fitRef = _jesFit->EvalPar(&ptref,par);
        //   // double fitRef = _jesFit->EvalPar(&pt,par);
        //   // cout << "in " << ptref << " " << pt << " " << data_point << " "
        //   << fitRef << " out " << data_point*fitRef << endl;
        //   // data_point *= fitRef;
        //   // sigma *= fitRef;
        // }

        // Calculate total shift caused by all nuisance parameters
        double shift = 0;
        for (auto [nuisance_name, nuisance] : nuisances) {
          if (dt_name != nuisance->appliesTo())
            continue;
          // if (FindInString("prob",nuisance_name.Data())
          // !=nuisance->appliesTo()) continue; if (debug) PrintLine("adding1
          // "+nuisance_name+" to shift for "+
          // dt_name+"("+nuisance->appliesTo()+")", yellow);
          TH1D *hnuis = nuisance->hist();
          assert(hnuis);
          shift += nuisance_pars[nuisance->index()] * hnuis->GetBinContent(hnuis->FindBin(pt));
        }

        // Add chi2 from residual
        // cout << "data_point: " << data_point << " shift: " << shift << " fit:
        // " << fit << " sigma: " << sigma << " err: " <<
        // oplus(sigma,GlobalFit::globalErrMin) << endl; if (debug)
        // PrintLine("chi2 before point:
        // "+to_string_with_precision(chi2),yellow);
        double fitPFG_delta = 0.25;
        double chi;
        if ("Resp" == dt_type) {
          chi = (data_point + shift - fit) / oplus(sigma, GlobalFit::globalErrMin);
        } else {
          chi = max(fabs(data_point - fit) - GlobalFit::ScaleFullSimShape * fitPFG_delta, 0.) / sigma;
          // cout << red << "FIND ME " << dt_type << " " << chi << " " << shift
          // << reset << endl;
        }

        chi2 += chi * chi;
        // cout << cyan << "Adding to chi2 " << chi * chi << " " << dt_name << "
        // " << pt << " " << data_point << " " << shift << " " << fit << " " <<
        // sigma << reset << endl; cout << cyan << "Adding to chi2 " <<
        // data_point << " " << sigma << " " << shift << " " << 0 << " " << chi
        // << reset << endl;
        ++nFittedDataPoints;
        //
        // // Store shifted data
        assert(graph_output->GetN() == graph_input->GetN() && graph_output->GetX()[bin] == pt);
        graph_output->SetPoint(bin, pt, ("Resp" == dt_type) ? data_point + shift : fit);
        graph_variation->SetPoint(bin, pt, ("Resp" == dt_type) ? shift : fit);

        // For multijets, store also downward extrapolation
        if (is_multijet_resp) {
          // TGraphErrors *recoil_raw       = dt_recoil->raw();
          // assert(recoil_raw);
          TGraphErrors *recoil_output = dt_recoil->output();
          assert(recoil_output);
          double ptref = recoil_output->GetX()[bin];
          double jesref = _jesFit->EvalPar(&ptref, par);
          double jes = _jesFit->EvalPar(&pt, par);
          // cout << "Shifting " << pt << " " << ptref << endl;
          recoil_output->SetPoint(bin, ptref, jes * jesref / (data_point + shift));
        } // multijet
      } // for point in graph
    } // for graph

    // Add chi2 from nuisance parameters
    // if (debug) PrintLine("chi2 before nuis:
    // "+to_string_with_precision(chi2),yellow);
    for (int ipar = 0; ipar < _nNuisancePars; ++ipar) {
      chi2 += nuisance_pars[ipar] * nuisance_pars[ipar];
      ++nFittedDataPoints;
    }
    // if (debug) PrintLine("chi2 before pars:
    // "+to_string_with_precision(chi2),yellow); Add penalty for fit parameters
    // (Bayesian prior, essentially)
    if (GlobalFit::penalizeFitPars) {
      for (int ipar = 0; ipar < _nFitPars; ++ipar) {
        chi2 += fit_pars[ipar] * fit_pars[ipar];
        ++nFittedDataPoints;
      }
    } // penalizeFitPars

    // if (debug) PrintLine("chi2 before next
    // point+to_string_with_precision(chi2),yellow);

    // Give some feedback on progress in case loop gets stuck
    if ((++fit_counter) % 100 == 0)
      cout << "." << flush;
  } // if flag
  else {
    if (grad) {
    }; // suppress warning;
    return;
  }

  // cout << yellow << "chi2: " << chi2 << reset << endl;

} // jesFitter
