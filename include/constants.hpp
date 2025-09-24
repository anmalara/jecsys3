#pragma once
#include "Containers.hpp"
#include "utils.hpp"

static const bool debug = true;
static constexpr double func_range_min = 10.;   // Define fitting range
static constexpr double func_range_max = 6500.; // Define fitting range
static constexpr double ptmin_multijet = -1;    // Min pt considered for multijet recoil
static constexpr double ptmax_multijet = 1300;  // Max pt considered for multijet recoil

static constexpr double ptmin_pf = 40;   // Min pt considered for pf composition
static constexpr double ptmax_pf = 1000; // Max pt considered for pf composition

static const std::map<TString, TString> input_fnames = {
    {"jes", "rootfiles/jecdataRUN.root"},
};

// clang-format off
static const std::map<TString, std::map<TString,TString>> reference_obj_map = {
  {"hjesref",       { {"fname", "jes"}, {"type", "jes"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_l2l3res"},}},
  {"herr",          { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr"},}},

  {"herr_ref",      { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_ref"},}},
  {"hrun1",         { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/hrun1"},}},
  {"hjes",          { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/hjes"},}},
  {"herr_l2l3res",  { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_l2l3res"},}},
  {"herr_ref",      { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_ref"},}},
  {"herr_noflv",    { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_noflv"},}},
  {"herr_spr",      { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_spr"},}},
  {"herr_pu",       { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_pu"},}},
};

static const std::map<TString, std::map<TString,TString>> input_hnames_map = {
  {"Resp_zjet_mpf",     { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/mpfchs1_zjet_a100"},}},
  {"Resp_zjet_db",      { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/ptchs_zjet_a100"},}},

  {"Resp_gamjet_mpf",   { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/mpfchs1_gamjet_a100"},}},
  {"Resp_gamjet_db",    { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/ptchs_gamjet_a100"},}},

  {"Resp_hadw_mpf",     { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/mpfchs1_hadw_a30"},}},
  // {"Resp_hadw_db",      { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/ptchs_hadw_a30"},}},
  // {"Resp_hadw_db",      { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/mpfchs1_hadw_a30"},}},

  {"Resp_multijet_mpf", { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/mpfchs1_multijet_a30"},}},
  {"Resp_multijet_db",  { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/ptchs_multijet_a30"},}},

  {"Resp_incljet_mpf",  { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/mpfchs1_incjet_a100"},}},
  {"Resp_incljet_db",   { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/ptchs_incjet_a100"},}},

  {"chf",               { {"fname", "jes"}, {"type", "chf"},  {"hname","MODE/etaETAMIN-ETAMAX/chf_cmb_ren"},}},
  {"nhf",               { {"fname", "jes"}, {"type", "nhf"},  {"hname","MODE/etaETAMIN-ETAMAX/nhf_cmb_ren"},}},
  {"nef",               { {"fname", "jes"}, {"type", "nef"},  {"hname","MODE/etaETAMIN-ETAMAX/nef_cmb_ren"},}},

};


static const std::map<TString, std::map<TString,TString>> recoil_hnames_map = {
  {"Resp_multijet_recoil_mpf", { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/crecoil_multijet_a30"},}},
  {"Resp_multijet_recoil_db",  { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/crecoil_multijet_a30"},}},
};

static const std::map<TString, std::map<TString,TString>> kfsr_hnames_map = {
  {"fsr_zjet_mpf",     { {"appliesTo", "Resp_zjet_mpf"},     {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_zjet"},}},
  {"fsr_gamjet_mpf",   { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_gamjet"},}},
  {"fsr_multijet_mpf", { {"appliesTo", "Resp_multijet_mpf"}, {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_multijet"},}},

  {"fsr_zjet_db",      { {"appliesTo", "Resp_zjet_db"},      {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_zjet"},}},
  {"fsr_gamjet_db",    { {"appliesTo", "Resp_gamjet_db"},    {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_gamjet"},}},
  {"fsr_multijet_db",  { {"appliesTo", "Resp_multijet_db"},  {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_multijet"},}},
};

static const std::map<TString, std::map<TString,TString>> nuisances_map = {
  {"Resp_uncl_zjet_mpf",        { {"appliesTo", "Resp_zjet_mpf"},     {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_zjet_mpfu1"},}},
  {"Resp_add_jet_zjet_mpf",     { {"appliesTo", "Resp_zjet_mpf"},     {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_zjet_mpfn1"},}},
  {"Resp_scale_zjet_mpf",       { {"appliesTo", "Resp_zjet_mpf"},     {"fname", "jes"}, {"hname","scale0.2"},}},
  {"Resp_hdmscale_zjet_mpf",    { {"appliesTo", "Resp_zjet_mpf"},     {"fname", "jes"}, {"hname","scale0.2"},}},

  // {"Resp_scale_zjet_db",       { {"appliesTo", "Resp_zjet_db"},     {"fname", "jes"}, {"hname","scale0.2"},}},
  // {"Resp_hdmscale_zjet_db",    { {"appliesTo", "Resp_zjet_db"},     {"fname", "jes"}, {"hname","scale0.2"},}},

  {"Resp_uncl_gamjet_mpf",      { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_gamjet_mpfu1"},}},
  {"Resp_add_jet_gamjet_mpf",   { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_gamjet_mpfn1"},}},
  {"Resp_scale_gamjet_mpf",     { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","scale0.5"},}},
  {"Resp_hdmscale_gamjet_mpf",  { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","scale0.2"},}},

  {"Resp_scale_gamjet_db",     { {"appliesTo", "Resp_gamjet_db"},   {"fname", "jes"}, {"hname","scale0.5"},}},
  {"Resp_hdmscale_gamjet_db",  { {"appliesTo", "Resp_gamjet_db"},   {"fname", "jes"}, {"hname","scale0.2"},}},

  {"Resp_zee_gamjet_eig0_mpf",  { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/sys/zee_gamjet_eig0"},}},
  {"Resp_zee_gamjet_eig1_mpf",  { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/sys/zee_gamjet_eig1"},}},
  {"Resp_zee_gamjet_eig2_mpf",  { {"appliesTo", "Resp_gamjet_mpf"},   {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/sys/zee_gamjet_eig2"},}},

  {"Resp_hadw_fitprob_mpf",     { {"appliesTo", "Resp_hadw_mpf"},     {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/sys/hadw_ptave_fitprob"},}},
  {"Resp_hadw_fitprob2_mpf",    { {"appliesTo", "Resp_hadw_mpf"},     {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/sys/hadw_ptave_fitprob2"},}},

  {"Resp_uncl_multijet_mpf",    { {"appliesTo", "Resp_multijet_mpf"}, {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_multijet_mpfu1"},}},
  {"Resp_add_jet_multijet_mpf", { {"appliesTo", "Resp_multijet_mpf"}, {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_multijet_mpfn1"},}},

  {"Resp_uncl_zjet_db",         { {"appliesTo", "Resp_zjet_db"},      {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_zjet_mpfu1"},}},
  {"Resp_add_jet_zjet_db",      { {"appliesTo", "Resp_zjet_db"},      {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_zjet_mpfn1"},}},
  {"Resp_uncl_gamjet_db",       { {"appliesTo", "Resp_gamjet_db"},    {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_gamjet_mpfu1"},}},
  {"Resp_add_jet_gamjet_db",    { {"appliesTo", "Resp_gamjet_db"},    {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_gamjet_mpfn1"},}},
  {"Resp_uncl_multijet_db",     { {"appliesTo", "Resp_multijet_db"},  {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_multijet_mpfu1"},}},
  {"Resp_add_jet_multijet_db",  { {"appliesTo", "Resp_multijet_db"},  {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_ptchs_multijet_mpfn1"},}},

};

static const std::map<TString, std::map<TString,TString>> shapes_map = {
    {"ecalm3",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalm3"}, {"appliesTo","Resp"}, {"form", "0.9734400837775052+log(x)*(0.01602168946078072+log(x)*(-0.0031069522045791215+log(x)*(-9.9662282335667e-05+log(x)*(3.461402158662177e-05+log(x)*(3.847832300063758e-06+log(x)*(-2.273285824254875e-07+log(x)*(-8.24313431758651e-08+log(x)*5.937419232146847e-09)))))))"},}},
    {"ecalm3_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalm3"}, {"appliesTo","chf"}, {"form", "-0.03123709391879015+log(x)*(0.03297298825067952+log(x)*(-0.010872334486666175+log(x)*(0.0010010717310248926+log(x)*(0.00015884108107291924+log(x)*(-2.2597079401917106e-05+log(x)*(-2.951789757144343e-06+log(x)*(5.894907780032335e-07+log(x)*-2.4883262238547824e-08)))))))"},}},
    {"ecalm3_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalm3"}, {"appliesTo","nhf"}, {"form", "0.01046207792614067+log(x)*(-0.009547425250785974+log(x)*(0.002880929360194274+log(x)*(-0.0001844551176742383+log(x)*(-5.4240462860280996e-05+log(x)*(9.413082549598995e-06+log(x)*-4.1121492972194564e-07)))))"},}},
    {"ecalm3_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalm3"}, {"appliesTo","nef"}, {"form", "0.02920526296840113+log(x)*(-0.03183662139386598+log(x)*(0.011014701620292305+log(x)*(-0.0012060488214486599+log(x)*(-0.00013067436403046382+log(x)*(2.5081953625741673e-05+log(x)*(2.4063714264936176e-06+log(x)*(-5.90319587531005e-07+log(x)*2.6919550498585114e-08)))))))"},}},

    {"ecalGain1p3",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain1p3"}, {"appliesTo","Resp"}, {"form", "0.9682414331248546+log(x)*(0.02675066236996259+log(x)*(-0.0062958907975048+log(x)*(-0.00011269731609060379+log(x)*(0.00012727927591589427+log(x)*(1.1824502595226483e-05+log(x)*(-1.4762697284780915e-06+log(x)*(-4.053098843629872e-07+log(x)*(-2.245802810203828e-08+log(x)*(4.389234213314463e-09+log(x)*(9.742292648743312e-10+log(x)*(5.2430414855565186e-11+log(x)*(-1.0759827288051267e-11+log(x)*(-2.1609864709101695e-12+log(x)*1.8968666228776667e-13)))))))))))))"},}},
    {"ecalGain1p3_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain1p3"}, {"appliesTo","chf"}, {"form", "0.010801781039921715+log(x)*(-0.009137263880398025+log(x)*(0.0020946378830883406+log(x)*(8.545058836619253e-05+log(x)*(-5.2912090893720865e-05+log(x)*(-4.45233778459218e-06+log(x)*(6.535592441378884e-07+log(x)*(1.5502886370466996e-07+log(x)*(7.540958023517135e-09+log(x)*(-1.7576835527074475e-09+log(x)*(-3.7868736100832013e-10+log(x)*(-2.0804958317598485e-11+log(x)*(4.348780928443551e-12+log(x)*(9.074722179897678e-13+log(x)*-8.072066052987796e-14)))))))))))))"},}},
    {"ecalGain1p3_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain1p3"}, {"appliesTo","nhf"}, {"form", "0.004611399364347392+log(x)*(-0.0037447771322590973+log(x)*(0.0008046225575780142+log(x)*(4.6022595174098984e-05+log(x)*(-2.1940067766724893e-05+log(x)*(-1.8471352524552112e-06+log(x)*(3.4428220731427714e-07+log(x)*(6.675663188164672e-08+log(x)*(-1.940149203435941e-09+log(x)*(-1.6646877835033754e-09+log(x)*1.1637751255658771e-10)))))))))"},}},
    {"ecalGain1p3_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain1p3"}, {"appliesTo","nef"}, {"form", "-0.009968502395309723+log(x)*(0.00831782899470734+log(x)*(-0.0018231402889229585+log(x)*(-0.00011077697433644632+log(x)*(5.045949058798476e-05+log(x)*(4.523585453944854e-06+log(x)*(-6.126746948533644e-07+log(x)*(-1.4620159510462396e-07+log(x)*(-6.919560432109076e-09+log(x)*(1.5994689854974068e-09+log(x)*(3.339698972039721e-10+log(x)*(1.7767546819317853e-11+log(x)*(-3.653225541199466e-12+log(x)*(-7.427548896616447e-13+log(x)*6.431088551925153e-14)))))))))))))"},}},

    {"ecalGain6p3",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain6p3"}, {"appliesTo","Resp"}, {"form", "0.7920242479718892+log(x)*(0.23322352164336912+log(x)*(-0.09081014150393805+log(x)*(0.011226261858872225+log(x)*(0.0010708561732756332+log(x)*(-0.00022122459654658086+log(x)*(-2.537606796800616e-05+log(x)*(2.051832665704375e-06+log(x)*(5.6111478232743e-07+log(x)*(2.0223401617940605e-08+log(x)*(-6.724149916716133e-09+log(x)*(-9.93570464340686e-10+log(x)*(2.2582102542591242e-11+log(x)*(2.0069138657934743e-11+log(x)*-1.2137615810730823e-12)))))))))))))"},}},
    {"ecalGain6p3_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain6p3"}, {"appliesTo","chf"}, {"form", "-0.016540270892427912+log(x)*(0.013747374858011842+log(x)*(-0.0029951136111960076+log(x)*(-0.0002027527047702164+log(x)*(9.410335723083026e-05+log(x)*(6.549110933150457e-06+log(x)*(-1.5171275753729308e-06+log(x)*(-2.4078399152625785e-07+log(x)*(1.0675182277024274e-08+log(x)*(5.752799533252303e-09+log(x)*-4.1672270783664126e-10)))))))))"},}},
    {"ecalGain6p3_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain6p3"}, {"appliesTo","nhf"}, {"form", "-0.003993119704412426+log(x)*(0.0038357595997896374+log(x)*(-0.0012077277709797592+log(x)*(8.025889334010765e-05+log(x)*(2.901922092564388e-05+log(x)*(-5.6295029505232884e-06+log(x)*2.8431406941244716e-07)))))"},}},
    {"ecalGain6p3_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain6p3"}, {"appliesTo","nef"}, {"form", "0.04185727246259771+log(x)*(-0.041520481591893144+log(x)*(0.013542951502934423+log(x)*(-0.0009760940134363541+log(x)*(-0.0002908808728757314+log(x)*(3.7187221593763936e-05+log(x)*(5.071731282281064e-06+log(x)*(-1.0658708453447075e-06+log(x)*4.7895714162402225e-08)))))))"},}},

    {"ecalGain12p3",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain12p3"}, {"appliesTo","Resp"}, {"form", "0.7425194340492929+log(x)*(0.2603168286775565+log(x)*(-0.08997803286397986+log(x)*(0.010292523072141471+log(x)*(0.000729107404861251+log(x)*(-0.0001743683892082836+log(x)*(-1.2080432123761862e-05+log(x)*(3.4593235450579746e-06+log(x)*-1.6061744156070108e-07)))))))"},}},
    {"ecalGain12p3_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain12p3"}, {"appliesTo","chf"}, {"form", "0.04908506444601239+log(x)*(-0.042800854254722936+log(x)*(0.010267569113912618+log(x)*(-2.4647415719618665e-05+log(x)*(-0.00019841476588200826+log(x)*(-9.70843216255186e-06+log(x)*(2.8906179768650086e-06+log(x)*(4.3937356548006303e-07+log(x)*(-1.6362949193849825e-08+log(x)*(-1.0257022100313807e-08+log(x)*6.993852298390174e-10)))))))))"},}},
    {"ecalGain12p3_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain12p3"}, {"appliesTo","nhf"}, {"form", "0.0010638937036792063+log(x)*(-0.0005839876837235162+log(x)*(8.454974536612714e-06+log(x)*(9.192372967334618e-06+log(x)*(1.4248145752591279e-07+log(x)*(-1.5933884231515725e-07+log(x)*9.482797531276756e-09)))))"},}},
    {"ecalGain12p3_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain12p3"}, {"appliesTo","nef"}, {"form", "-0.05221750600313704+log(x)*(0.046025993827596666+log(x)*(-0.011472364820063873+log(x)*(0.00021735149031875008+log(x)*(0.00020179592854752167+log(x)*(6.306690242776466e-06+log(x)*(-3.001164052611249e-06+log(x)*(-3.8336232684914437e-07+log(x)*(2.0024423679326965e-08+log(x)*(9.121078625540033e-09+log(x)*-6.467888132111031e-10)))))))))"},}},

    {"hadHcalp3",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalp3"}, {"appliesTo","Resp"}, {"form", "0.9453492957483248+log(x)*(0.06934541092828847+log(x)*(-0.03379584216240194+log(x)*(0.008505431802514006+log(x)*(-0.0011360503817276356+log(x)*(7.724380629442239e-05+log(x)*-2.113537847670741e-06)))))"},}},
    {"hadHcalp3_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalp3"}, {"appliesTo","chf"}, {"form", "0.022268438721858727+log(x)*(-0.022203171571895847+log(x)*(0.007074480059271589+log(x)*(-0.000641969821834951+log(x)*(-0.00010939676058479088+log(x)*(2.36321626515414e-05+log(x)*-1.1453261916508196e-06)))))"},}},
    {"hadHcalp3_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalp3"}, {"appliesTo","nhf"}, {"form", "0.04630781898146411+log(x)*(-0.03771272333197105+log(x)*(0.009121984095758411+log(x)*(3.3991124726973026e-05+log(x)*(-0.00017109691027934764+log(x)*(-9.474838209638437e-06+log(x)*(2.3590505039502995e-06+log(x)*(3.7475935201232643e-07+log(x)*(-1.2146782911117676e-08+log(x)*(-8.346877459548043e-09+log(x)*5.53895729791253e-10)))))))))"},}},
    {"hadHcalp3_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalp3"}, {"appliesTo","nef"}, {"form", "-0.031718556668351966+log(x)*(0.02585307309199026+log(x)*(-0.00583788861024489+log(x)*(-0.00013181215715265138+log(x)*(0.00010773139401963508+log(x)*(9.654149004914184e-06+log(x)*(-8.833241689301541e-07+log(x)*(-2.592928006500091e-07+log(x)*(-1.5997357718612633e-08+log(x)*(2.6392395059534467e-09+log(x)*(5.65964189452409e-10+log(x)*(-4.819817043122091e-11)))))))))))"},}},

    {"hadHcalH097",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalH097"}, {"appliesTo","Resp"}, {"form", "1.2064948461463798+log(x)*(-0.21339276784465533+log(x)*(0.07362452719452818+log(x)*(-0.006861126865705539+log(x)*(-0.001031746519025494+log(x)*(0.00011452064645044169+log(x)*(2.1722659432461185e-05+log(x)*(-6.315388919356715e-07+log(x)*(-3.9559138135721247e-07+log(x)*(-1.162968239527087e-08+log(x)*(7.026507945386855e-09+log(x)*(-3.423667621251371e-10)))))))))))"},}},
    {"hadHcalH097_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalH097"}, {"appliesTo","chf"}, {"form", "-0.09841011745023989+log(x)*(0.09729220605997899+log(x)*(-0.030646359323668184+log(x)*(0.001869349270121861+log(x)*(0.0005544358921908957+log(x)*(-2.062337681864308e-05+log(x)*(-1.096480059926944e-05+log(x)*(-5.41035636928893e-07+log(x)*(1.3497133028572785e-07+log(x)*(2.214600730825385e-08+log(x)*(-4.4685399318640276e-10+log(x)*(-4.61189358564845e-10+log(x)*2.826042964152188e-11)))))))))))"},}},
    {"hadHcalH097_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalH097"}, {"appliesTo","nhf"}, {"form", "0.020560328659924337+log(x)*(-0.020339445725213175+log(x)*(0.005852159304843714+log(x)*(-9.681067673068904e-05+log(x)*(-0.0001341805984587659+log(x)*(-3.886789905486341e-06+log(x)*(2.4409530952770295e-06+log(x)*(2.965384033333383e-07+log(x)*(-2.1691271880242577e-08+log(x)*(-8.018022930635208e-09+log(x)*6.275428667093182e-10)))))))))"},}},
    {"hadHcalH097_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalH097"}, {"appliesTo","nef"}, {"form", "0.054135418076619105+log(x)*(-0.05447703648087184+log(x)*(0.01769759116910445+log(x)*(-0.0011538638113039287+log(x)*(-0.00034748636580395683+log(x)*(1.7966486144085983e-05+log(x)*(7.249082940812804e-06+log(x)*(2.47532189435327e-07+log(x)*(-1.0266361907599065e-07+log(x)*(-1.4251430532367734e-08+log(x)*(5.641310006712985e-10+log(x)*(3.245853548852167e-10+log(x)*-2.189194551912771e-11)))))))))))"},}},
    
    {"hadHcalH100",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalH100"}, {"appliesTo","Resp"}, {"form", "1.130052661789372+log(x)*(-0.13079147907799865+log(x)*(0.04276585917989275+log(x)*(-0.003332436950773344+log(x)*(-0.0006187857662424555+log(x)*(3.79400649113919e-05+log(x)*(1.1905556453697553e-05+log(x)*(3.93014364954051e-07+log(x)*(-1.4464411312124024e-07+log(x)*(-1.996902481527348e-08+log(x)*(5.559456115522524e-10+log(x)*(3.9330884826947954e-10+log(x)*-2.3690755912284858e-11)))))))))))"},}},
    {"hadHcalH100_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalH100"}, {"appliesTo","chf"}, {"form", "-0.13461891719050267+log(x)*(0.13950747209325315+log(x)*(-0.04838986792645597+log(x)*(0.004573546528807363+log(x)*(0.0006780710752729117+log(x)*(-7.740182453018314e-05+log(x)*(-1.4539101313980487e-05+log(x)*(4.498236014478943e-07+log(x)*(2.7331891377833735e-07+log(x)*(7.844390546665349e-09+log(x)*(-5.038367881391573e-09+log(x)*(2.5251790284821505e-10)))))))))))"},}},
    {"hadHcalH100_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalH100"}, {"appliesTo","nef"}, {"form", "0.02309247999552963+log(x)*(-0.022220913070106783+log(x)*(0.006281059226005439+log(x)*(-0.0001280069677541228+log(x)*(-0.0001332459363043535+log(x)*(-3.2524998492570707e-06+log(x)*(2.3071844515312935e-06+log(x)*(2.675698195383213e-07+log(x)*(-1.954784438527377e-08+log(x)*(-6.994423593042257e-09+log(x)*5.37713782107747e-10)))))))))"},}},
    {"hadHcalH100_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalH100"}, {"appliesTo","nhf"}, {"form", "0.08605941248699643+log(x)*(-0.09187959536761778+log(x)*(0.03345316283474071+log(x)*(-0.003550492925799591+log(x)*(-0.00045812412858673107+log(x)*(6.834235204510362e-05+log(x)*(1.0142782963010441e-05+log(x)*(-5.979345998485299e-07+log(x)*(-2.1148223286798495e-07+log(x)*(-2.3437847060678675e-09+log(x)*(4.045894692359442e-09+log(x)*(-2.2426213552502653e-10)))))))))))"},}},

    {"hadHcalH106",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalH106"}, {"appliesTo","Resp"}, {"form", "1.2002396997187292+log(x)*(-0.208534508824019+log(x)*(0.07467793253670006+log(x)*(-0.008399866862984912+log(x)*(-0.0007916776288672237+log(x)*(0.00017058298855723668+log(x)*(1.3474072594322465e-05+log(x)*(-3.709358352527495e-06+log(x)*1.7406693572463496e-07)))))))"},}},
    {"hadHcalH106_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalH106"}, {"appliesTo","chf"}, {"form", "-0.08433687836218724+log(x)*(0.08109487222754519+log(x)*(-0.02445447086648641+log(x)*(0.0012593702395378732+log(x)*(0.0004417014018475422+log(x)*(-1.0040800089366495e-05+log(x)*(-8.073579557808632e-06+log(x)*(-4.958700686070162e-07+log(x)*(8.450897079982737e-08+log(x)*(1.5913012042505518e-08+log(x)*(-8.789432752850473e-11+log(x)*(-3.001817261197809e-10+log(x)*1.6609088119567825e-11)))))))))))"},}},
    {"hadHcalH106_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalH106"}, {"appliesTo","nef"}, {"form", "0.01682050661553866+log(x)*(-0.013817580036767587+log(x)*(0.0028173989211748417+log(x)*(0.00018456131214613768+log(x)*(-5.3615695038616946e-05+log(x)*(-7.871751410516311e-06+log(x)*(1.2384847637067754e-07+log(x)*(1.5510749889006076e-07+log(x)*(1.9503736482154218e-08+log(x)*(1.3351121611001442e-10+log(x)*(-3.263803742916996e-10+log(x)*(-4.40170047030706e-11+log(x)*4.960294021819479e-12)))))))))))"},}},
    {"hadHcalH106_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalH106"}, {"appliesTo","nhf"}, {"form", "0.07177245586718797+log(x)*(-0.07054860287331657+log(x)*(0.022055221539709513+log(x)*(-0.0011975824129994977+log(x)*(-0.0004587848089689023+log(x)*(1.8079856513327183e-05+log(x)*(9.149001589746603e-06+log(x)*(3.773331556312673e-07+log(x)*(-1.1989375087837067e-07+log(x)*(-1.7610451041721294e-08+log(x)*(5.438532856228997e-10+log(x)*(3.7438082673291687e-10+log(x)*-2.389511810683862e-11)))))))))))"},}},
    
    {"trkEffNtrk1m3",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk1m3"}, {"appliesTo","Resp"}, {"form", "0.9705542049351432+log(x)*(0.00833644189700795+log(x)*(-0.0005268782538820215+log(x)*(-5.830565847569534e-06+log(x)*(1.3699865668494707e-06+log(x)*(-1.6432688595714373e-07+log(x)*(-5.0929660282295565e-08+log(x)*(-4.523210064030283e-09+log(x)*(1.3538693431885973e-10+log(x)*(1.0028986714328572e-10+log(x)*(1.6207126119192183e-11+log(x)*(1.3982431483482965e-12+log(x)*(-1.0410078239872199e-15+log(x)*(-2.331308768600071e-14+log(x)*(-4.240671896800829e-15+log(x)*(-2.8921802625829143e-16+log(x)*6.593562793879253e-17)))))))))))))))"},}},
    {"trkEffNtrk1m3_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk1m3"}, {"appliesTo","chf"}, {"form", "0.13422523845324447+log(x)*(-0.10529144710632915+log(x)*(0.021285453229086246+log(x)*(0.0005303470903335295+log(x)*(-0.0003205488321033704+log(x)*(-3.0189773329398807e-05+log(x)*(2.151591048865962e-06+log(x)*(7.086015364856973e-07+log(x)*(5.3126601004137086e-08+log(x)*(-4.299785918200419e-09+log(x)*(-1.4172317282916526e-09+log(x)*(-1.1232020106648615e-10+log(x)*(1.1338515305363244e-11+log(x)*(3.2961385488806393e-12+log(x)*-2.484050304888955e-13)))))))))))))"},}},
    {"trkEffNtrk1m3_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk1m3"}, {"appliesTo","nef"}, {"form", "-0.06784427539441518+log(x)*(0.05211631720565596+log(x)*(-0.010304588456854735+log(x)*(-0.00029870373257875644+log(x)*(0.0001541599726571226+log(x)*(1.5460122296101724e-05+log(x)*(-9.494645720833492e-07+log(x)*(-3.484603934279158e-07+log(x)*(-2.7848336476245492e-08+log(x)*(1.9261271214685895e-09+log(x)*(7.027221346266921e-10+log(x)*(5.865670800340777e-11+log(x)*(-5.389347321131359e-12+log(x)*(-1.6695623137951887e-12+log(x)*1.235754836701425e-13)))))))))))))"},}},
    {"trkEffNtrk1m3_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk1m3"}, {"appliesTo","nhf"}, {"form", "-0.06525627279159024+log(x)*(0.05257251106216328+log(x)*(-0.011005092112109021+log(x)*(-0.00018038851414241654+log(x)*(0.00016362960634284284+log(x)*(1.3456450002031199e-05+log(x)*(-1.2297249382739693e-06+log(x)*(-3.3698635878931567e-07+log(x)*(-2.211687738234804e-08+log(x)*(2.317822923311861e-09+log(x)*(6.46656690285443e-10+log(x)*(4.6176727385714804e-11+log(x)*(-5.442755823246361e-12+log(x)*(-1.4140270047410726e-12+log(x)*1.0906746624434398e-13)))))))))))))"},}},

    {"trkEffNtrk2ToInfm3",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk2ToInfm3"}, {"appliesTo","Resp"}, {"form", "1.0093778477781266+log(x)*(-0.0035650903137754517+log(x)*(-0.0006013033855872676+log(x)*(0.0001400714543912632+log(x)*(1.1161943432170985e-05+log(x)*(-9.63349516736549e-07+log(x)*(-1.866580573001355e-07+log(x)*(5.877876234170427e-10+log(x)*1.145614695638029e-09)))))))"},}},
    {"trkEffNtrk2ToInfm3_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk2ToInfm3"}, {"appliesTo","chf"}, {"form", "-0.04454858118799183+log(x)*(0.03231264555268816+log(x)*(-0.005899837296451476+log(x)*(-0.0003712999064835149+log(x)*(8.081230684601325e-05+log(x)*(1.0305144994028992e-05+log(x)*(-5.192327989378855e-07+log(x)*(-2.0209884080282563e-07+log(x)*1.3989897284832443e-08)))))))"},}},
    {"trkEffNtrk2ToInfm3_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk2ToInfm3"}, {"appliesTo","nef"}, {"form", "0.1744087779610073+log(x)*(-0.16554856719677616+log(x)*(0.05300884958884678+log(x)*(-0.0047974387482971085+log(x)*(-0.0006274059048429867+log(x)*(9.689895087842393e-05+log(x)*(1.0129400189717513e-05+log(x)*(-2.2281211734858754e-06+log(x)*9.738599611522415e-08)))))))"},}},
    {"trkEffNtrk2ToInfm3_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk2ToInfm3"}, {"appliesTo","nhf"}, {"form", "-0.09471451642969325+log(x)*(0.1018306637543447+log(x)*(-0.03840950932345122+log(x)*(0.004961106932913725+log(x)*(0.000290719242281886+log(x)*(-8.534770951967594e-05+log(x)*(-5.1451312648149e-06+log(x)*(1.6535142121500972e-06+log(x)*-7.826902055016421e-08)))))))"},}},

    {"trkEff0999Nm1",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0999Nm1"}, {"appliesTo","Resp"}, {"form", "1.0013242599898862+log(x)*(-0.0010795263994483448+log(x)*(0.00011741219229311441+log(x)*(4.0749242046975554e-05+log(x)*(-8.76206874199126e-06+log(x)*(-1.0688397421919224e-06+log(x)*(8.995689670715133e-08+log(x)*(2.981216398525571e-08+log(x)*(1.876435727134481e-09+log(x)*(-3.37687091456301e-10+log(x)*(-7.250769091103905e-11+log(x)*(6.437052730060625e-12)))))))))))"},}},
    {"trkEff0999Nm1_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0999Nm1"}, {"appliesTo","chf"}, {"form", "0.009872286117708575+log(x)*(-0.009007607776654783+log(x)*(0.002242799061775209+log(x)*(2.326626380834588e-05+log(x)*(-4.805428815411052e-05+log(x)*(-3.8645742472605606e-06+log(x)*(4.1516149399019526e-07+log(x)*(1.0904730512630657e-07+log(x)*(6.334414408739287e-09+log(x)*(-1.0647853963868285e-09+log(x)*(-2.191915009203179e-10+log(x)*(1.8284440520896096e-11)))))))))))"},}},
    {"trkEff0999Nm1_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0999Nm1"}, {"appliesTo","nef"}, {"form", "0.006172731082289809+log(x)*(-0.0036624862755707194+log(x)*(0.0004126722221112118+log(x)*(7.543854311739134e-05+log(x)*(-2.43427660982336e-06+log(x)*(-1.4531594242577402e-06+log(x)*(-1.3638895607190525e-07+log(x)*(4.2564320071212e-09+log(x)*(2.8413253746746024e-09+log(x)*(3.9205327428660363e-10+log(x)*(1.9444682305515457e-11+log(x)*(-3.2162641192845514e-12+log(x)*(-8.739283031824556e-13+log(x)*(-7.558836209143028e-14+log(x)*1.2835643778023857e-14)))))))))))))"},}},
    {"trkEff0999Nm1_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0999Nm1"}, {"appliesTo","nhf"}, {"form", "-0.009216817682272335+log(x)*(0.007137480763116881+log(x)*(-0.0014133973086683355+log(x)*(-6.474776709052071e-05+log(x)*(2.259415847984751e-05+log(x)*(3.330108940806109e-06+log(x)*(8.864911751543125e-08+log(x)*(-3.6892918003283806e-08+log(x)*(-7.2542231906868926e-09+log(x)*(-6.69922019366393e-10+log(x)*(-6.300174599731578e-12+log(x)*(9.234681718888957e-12+log(x)*(1.6928571308091924e-12+log(x)*(1.1283873148815965e-13+log(x)*-2.4934664866337097e-14)))))))))))))"},}},

  {"const",    {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","const"},   {"appliesTo","Resp"}, {"form", "1"},}},
  {"ftd",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ftd"},     {"appliesTo","Resp"}, {"form", "-0.116-0.6417*pow(x/208.,-0.3051)+23.63/x"},}},
  {"fp",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fp"},      {"appliesTo","Resp"}, {"form", "-0.8295"},}},
  {"fhx",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhx"},     {"appliesTo","Resp"}, {"form", "0.8904+1.082*pow(x/1408,1.204)/(1+pow(x/1408,1.204))*(1-pow(x/1408,-1.204))"},}},
  {"fhh",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhh"},     {"appliesTo","Resp"}, {"form", "-0.7938-0.5798*pow(x/396.1,1.412)/(1+pow(x/396.1,1.412))*(1-pow(x/396.1,-1.412))"},}},
//   {"feh",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","feh"},     {"appliesTo","Resp"}, {"form", "-0.2603-0.2196*pow(x/409.4,1.276)/(1+pow(x/409.4,1.276))*(1-pow(x/409.4,-1.276))"},}},
  {"fhw",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhw"},     {"appliesTo","Resp"}, {"form", "0.3*(0.9526-0.3883*(1+(pow(x/1285,2.46)-1)/(pow(x/1285,2.46)+1))+18.1/x-2.062*log(x)/x)"},}},
  {"fl1",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fl1"},     {"appliesTo","Resp"}, {"form", "(1-(0.350077+0.553560*log(x)-0.0527681*pow(log(x),2))/x-1)"},}},
  {"ftd-ftm",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ftd-ftm"}, {"appliesTo","Resp"}, {"form", "3*((-0.116-0.6417*pow(x/208.,-0.3051)+23.63/x)-(0.2683-0.6994*pow(x/208.,-0.3051)+18.49/x))"},}},
  {"f1q3-1",   {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","f1q3-1"},  {"appliesTo","Resp"}, {"form", "0.01*(0.7966+0.9311*(pow(0.01*x,-1)-1))"},}},

  {"ftd_chf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ftd"}, {"appliesTo","chf"}, {"form", "1.982-2.678*(1+(pow(x/47.02,0.262)-1)/(pow(x/47.02,0.262)+1))+0.1494*pow(x,+0.3)-3.097/x"},}},
  {"ftd_nhf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ftd"}, {"appliesTo","nhf"}, {"form", "-0.01022-0.1962*(1+(pow(x/4000,3.071)-1)/(pow(x/4000,3.071)+1))+0.04211*pow(x,+0.3)+0.01005/x"},}},
  {"ftd_nef",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ftd"}, {"appliesTo","nef"}, {"form", "0.07453+0.1457*(1+(pow(x/1131,-3.68)-1)/(pow(x/1131,-3.68)+1))-0.4155*pow(x,-0.3)-1.878/x"},}},

  {"fp_chf",   {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fp"}, {"appliesTo","chf"}, {"form", "0.3333+0.7433*(1+(pow(x/1023,0.3926)-1)/(pow(x/1023,0.3926)+1))-0.09446*pow(x,0.2883)"},}},
  {"fp_nhf",   {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fp"}, {"appliesTo","nhf"}, {"form", "0.07395+0*(1+(pow(x/1000,258.2)-1)/(pow(x/1000,258.2)+1))+1.223e-05*pow(x,1.158)"},}},
  {"fp_nef",   {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fp"}, {"appliesTo","nef"}, {"form", "2.283+0*(1+(pow(x/1000,1.302)-1)/(pow(x/1000,1.302)+1))-2.738*pow(x,0.002452)"},}},

  {"fhx_chf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhx"}, {"appliesTo","chf"}, {"form", "-0.0637-0.2811*(1+(pow(x/4531,-0.3172)-1)/(pow(x/4531,-0.3172)+1))+1.071*pow(x,-0.153)"},}},
  {"fhx_nhf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhx"}, {"appliesTo","nhf"}, {"form", "-0.295+0.09444*(1+(pow(x/2713,0.06437)-1)/(pow(x/2713,0.06437)+1))+[p4]*pow(x,0.2845)"},}},
  {"fhx_nef",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhx"}, {"appliesTo","nef"}, {"form", "0.05474-0.003141*(1+(pow(x/798.6,78.84)-1)/(pow(x/798.6,78.84)+1))-0.000957*pow(x,0.76)"},}},

  {"fhh_chf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhh"}, {"appliesTo","chf"}, {"form", "0.1552-0.04221*(1+(pow(x/315,2.787)-1)/(pow(x/315,2.787)+1))-0.06628*pow(x,-0.2572)"},}},
  {"fhh_nhf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhh"}, {"appliesTo","nhf"}, {"form", "-0.2746-0.6358*(1+(pow(x/9664,0.6547)-1)/(pow(x/9664,0.6547)+1))+0.05559*pow(x,0.1816)"},}},
  {"fhh_nef",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhh"}, {"appliesTo","nef"}, {"form", "0.4158-2.14*(1+(pow(x/9426,0.1723)-1)/(pow(x/9426,0.1723)+1))+0.4111*pow(x,0.1937)"},}},

// //   {"feh_chf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","feh"}, {"appliesTo","chf"}, {"form", "0.06085+0*(1+(pow(x/1000,1.3)-1)/(pow(x/1000,1.3)+1))-0.008137*pow(x,0.2135)"},}},
// //   {"feh_nhf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","feh"}, {"appliesTo","nhf"}, {"form", "-0.03458+0*(1+(pow(x/1713,274.8)-1)/(pow(x/1713,274.8)+1))+0.01665*pow(x,0.2426)"},}},
// //   {"feh_nef",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","feh"}, {"appliesTo","nef"}, {"form", "-0.02364+0*(1+(pow(x/1481,246.2)-1)/(pow(x/1481,246.2)+1))-0.009737*pow(x,0.2576)"},}},

  {"fhw_chf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhw"}, {"appliesTo","chf"}, {"form", "-0.2176+1.064e-05*pow(x,1.373)+0/x"},}},
  {"fhw_nhf",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhw"}, {"appliesTo","nhf"}, {"form", "-5.151+4.495*pow(x,0.03335)-12.3/x"},}},
  {"fhw_nef",  {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","fhw"}, {"appliesTo","nef"}, {"form", "0.8417-0.2605*pow(x,0.2289)+2.426/x"},}},


  {"hadHcalp3",              {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalp3"},          {"appliesTo","Resp"}, {"form", "+2.792e+00-2.377e+00*log(x)+6.981e-01*pow(log(x),2)-7.477e-02*pow(log(x),3)+2.779e-03*pow(log(x),4)"},}},
  {"hadHcalp3_chf",          {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalp3"},          {"appliesTo","chf"},  {"form", "-2.368e+00+1.749e+00*log(x)-4.589e-01*pow(log(x),2)+4.307e-02*pow(log(x),3)-1.183e-03*pow(log(x),4)"},}},
  {"hadHcalp3_nhf",          {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalp3"},          {"appliesTo","nhf"},  {"form", "-4.683e+00+3.771e+00*log(x)-1.038e+00*pow(log(x),2)+1.224e-01*pow(log(x),3)-5.071e-03*pow(log(x),4)"},}},
  {"hadHcalp3_nef",          {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalp3"},          {"appliesTo","nef"},  {"form", "+8.141e+00-6.447e+00*log(x)+1.779e+00*pow(log(x),2)-2.020e-01*pow(log(x),3)+7.952e-03*pow(log(x),4)"},}},

  {"hadHcalZB097",           {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB097"},       {"appliesTo","Resp"}, {"form", "-7.714e+00+7.581e+00*log(x)-2.665e+00*pow(log(x),2)+3.940e-01*pow(log(x),3)-2.049e-02*pow(log(x),4)"},}},
  {"hadHcalZB097_chf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB097"},       {"appliesTo","chf"},  {"form", "+3.089e+00-2.918e+00*log(x)+1.006e+00*pow(log(x),2)-1.477e-01*pow(log(x),3)+7.648e-03*pow(log(x),4)"},}},
  {"hadHcalZB097_nhf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB097"},       {"appliesTo","nhf"},  {"form", "-9.328e+00+7.999e+00*log(x)-2.493e+00*pow(log(x),2)+3.336e-01*pow(log(x),3)-1.605e-02*pow(log(x),4)"},}},
  {"hadHcalZB097_nef",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB097"},       {"appliesTo","nef"},  {"form", "+4.981e+00-4.021e+00*log(x)+1.166e+00*pow(log(x),2)-1.446e-01*pow(log(x),3)+6.489e-03*pow(log(x),4)"},}},

  {"hadHcalZB100",           {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB100"},       {"appliesTo","Resp"}, {"form", "-9.705e+00+8.983e+00*log(x)-2.986e+00*pow(log(x),2)+4.190e-01*pow(log(x),3)-2.066e-02*pow(log(x),4)"},}},
  {"hadHcalZB100_chf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB100"},       {"appliesTo","chf"},  {"form", "+4.857e+00-4.366e+00*log(x)+1.427e+00*pow(log(x),2)-1.990e-01*pow(log(x),3)+9.812e-03*pow(log(x),4)"},}},
  {"hadHcalZB100_nhf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB100"},       {"appliesTo","nhf"},  {"form", "-9.351e+00+7.887e+00*log(x)-2.408e+00*pow(log(x),2)+3.141e-01*pow(log(x),3)-1.461e-02*pow(log(x),4)"},}},
  {"hadHcalZB100_nef",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB100"},       {"appliesTo","nef"},  {"form", "+3.789e+00-2.934e+00*log(x)+8.056e-01*pow(log(x),2)-9.292e-02*pow(log(x),3)+3.776e-03*pow(log(x),4)"},}},

  {"hadHcalZB106",           {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB106"},       {"appliesTo","Resp"}, {"form", "-2.023e+01+1.720e+01*log(x)-5.226e+00*pow(log(x),2)+6.681e-01*pow(log(x),3)-2.976e-02*pow(log(x),4)"},}},
  {"hadHcalZB106_chf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB106"},       {"appliesTo","chf"},  {"form", "+1.035e+01-8.927e+00*log(x)+2.779e+00*pow(log(x),2)-3.677e-01*pow(log(x),3)+1.719e-02*pow(log(x),4)"},}},
  {"hadHcalZB106_nhf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB106"},       {"appliesTo","nhf"},  {"form", "-1.206e+01+9.880e+00*log(x)-2.898e+00*pow(log(x),2)+3.578e-01*pow(log(x),3)-1.535e-02*pow(log(x),4)"},}},
  {"hadHcalZB106_nef",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","hadHcalZB106"},       {"appliesTo","nef"},  {"form", "+3.033e+00-2.060e+00*log(x)+4.513e-01*pow(log(x),2)-3.274e-02*pow(log(x),3)+1.314e-04*pow(log(x),4)"},}},

  {"ecalm3",                 {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalm3"},             {"appliesTo","Resp"}, {"form", "-9.628e-01+1.009e+00*log(x)-4.167e-01*pow(log(x),2)+5.326e-02*pow(log(x),3)-2.206e-03*pow(log(x),4)"},}},
  {"ecalm3_chf",             {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalm3"},             {"appliesTo","chf"},  {"form", "+5.813e+00-4.086e+00*log(x)+1.102e+00*pow(log(x),2)-1.220e-01*pow(log(x),3)+4.691e-03*pow(log(x),4)"},}},
  {"ecalm3_nhf",             {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalm3"},             {"appliesTo","nhf"},  {"form", "-1.636e-01+1.183e-01*log(x)-9.015e-03*pow(log(x),2)-3.084e-03*pow(log(x),3)+3.948e-04*pow(log(x),4)"},}},
  {"ecalm3_nef",             {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalm3"},             {"appliesTo","nef"},  {"form", "-5.823e+00+4.103e+00*log(x)-1.133e+00*pow(log(x),2)+1.300e-01*pow(log(x),3)-5.306e-03*pow(log(x),4)"},}},

  {"ecalGain1p3",            {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain1p3"},        {"appliesTo","Resp"}, {"form", "-8.351e+00+6.497e+00*log(x)-1.773e+00*pow(log(x),2)+1.998e-01*pow(log(x),3)-7.719e-03*pow(log(x),4)"},}},
  {"ecalGain1p3_chf",        {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain1p3"},        {"appliesTo","chf"},  {"form", "+2.612e+00-2.118e+00*log(x)+6.084e-01*pow(log(x),2)-7.302e-02*pow(log(x),3)+3.065e-03*pow(log(x),4)"},}},
  {"ecalGain1p3_nhf",        {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain1p3"},        {"appliesTo","nhf"},  {"form", "-9.112e-01+8.269e-01*log(x)-2.729e-01*pow(log(x),2)+3.890e-02*pow(log(x),3)-2.026e-03*pow(log(x),4)"},}},
  {"ecalGain1p3_nef",        {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain1p3"},        {"appliesTo","nef"},  {"form", "-3.583e-01+2.078e-01*log(x)-2.141e-02*pow(log(x),2)-4.756e-03*pow(log(x),3)+7.013e-04*pow(log(x),4)"},}},

  {"ecalGain6p3",            {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain6p3"},        {"appliesTo","Resp"}, {"form", "-5.587e+00+4.992e+00*log(x)-1.596e+00*pow(log(x),2)+2.157e-01*pow(log(x),3)-1.033e-02*pow(log(x),4)"},}},
  {"ecalGain6p3_chf",        {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain6p3"},        {"appliesTo","chf"},  {"form", "+2.414e+00-2.133e+00*log(x)+6.755e-01*pow(log(x),2)-9.053e-02*pow(log(x),3)+4.308e-03*pow(log(x),4)"},}},
  {"ecalGain6p3_nhf",        {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain6p3"},        {"appliesTo","nhf"},  {"form", "+7.040e-01-5.838e-01*log(x)+1.728e-01*pow(log(x),2)-2.155e-02*pow(log(x),3)+9.501e-04*pow(log(x),4)"},}},
  {"ecalGain6p3_nef",        {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain6p3"},        {"appliesTo","nef"},  {"form", "-3.593e+00+3.112e+00*log(x)-9.667e-01*pow(log(x),2)+1.272e-01*pow(log(x),3)-5.959e-03*pow(log(x),4)"},}},

  {"ecalGain12p3",           {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain12p3"},       {"appliesTo","Resp"}, {"form", "+1.542e+01-1.289e+01*log(x)+3.895e+00*pow(log(x),2)-4.819e-01*pow(log(x),3)+2.084e-02*pow(log(x),4)"},}},
  {"ecalGain12p3_chf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain12p3"},       {"appliesTo","chf"},  {"form", "-1.059e+01+8.124e+00*log(x)-2.321e+00*pow(log(x),2)+2.773e-01*pow(log(x),3)-1.169e-02*pow(log(x),4)"},}},
  {"ecalGain12p3_nhf",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain12p3"},       {"appliesTo","nhf"},  {"form", "-9.448e-02+1.776e-02*log(x)-1.417e-03*pow(log(x),2)-4.896e-04*pow(log(x),3)+5.913e-05*pow(log(x),4)"},}},
  {"ecalGain12p3_nef",       {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","ecalGain12p3"},       {"appliesTo","nef"},  {"form", "+1.063e+01-8.107e+00*log(x)+2.314e+00*pow(log(x),2)-2.758e-01*pow(log(x),3)+1.158e-02*pow(log(x),4)"},}},

  {"trkEff0999Nm1",          {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0999Nm1"},      {"appliesTo","Resp"}, {"form", "-2.415e+00+2.029e+00*log(x)-6.003e-01*pow(log(x),2)+7.296e-02*pow(log(x),3)-3.117e-03*pow(log(x),4)"},}},
  {"trkEff0999Nm1_chf",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0999Nm1"},      {"appliesTo","chf"},  {"form", "-1.755e+00+9.661e-01*log(x)-8.482e-02*pow(log(x),2)-2.171e-02*pow(log(x),3)+2.382e-03*pow(log(x),4)"},}},
  {"trkEff0999Nm1_nhf",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0999Nm1"},      {"appliesTo","nhf"},  {"form", "-2.462e+00+2.152e+00*log(x)-6.830e-01*pow(log(x),2)+9.276e-02*pow(log(x),3)-4.384e-03*pow(log(x),4)"},}},
  {"trkEff0999Nm1_nef",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0999Nm1"},      {"appliesTo","nef"},  {"form", "+5.553e+00-4.180e+00*log(x)+1.069e+00*pow(log(x),2)-1.073e-01*pow(log(x),3)+3.559e-03*pow(log(x),4)"},}},

  {"trkEff0998Nm1",          {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0998Nm1"},      {"appliesTo","Resp"}, {"form", "-4.645e+00+3.917e+00*log(x)-1.161e+00*pow(log(x),2)+1.410e-01*pow(log(x),3)-6.012e-03*pow(log(x),4)"},}},
  {"trkEff0998Nm1_chf",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0998Nm1"},      {"appliesTo","chf"},  {"form", "-5.411e+00+3.418e+00*log(x)-5.853e-01*pow(log(x),2)+6.119e-03*pow(log(x),3)+2.638e-03*pow(log(x),4)"},}},
  {"trkEff0998Nm1_nhf",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0998Nm1"},      {"appliesTo","nhf"},  {"form", "-5.007e+00+4.371e+00*log(x)-1.388e+00*pow(log(x),2)+1.887e-01*pow(log(x),3)-8.947e-03*pow(log(x),4)"},}},
  {"trkEff0998Nm1_nef",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEff0998Nm1"},      {"appliesTo","nef"},  {"form", "+1.121e+01-8.435e+00*log(x)+2.159e+00*pow(log(x),2)-2.174e-01*pow(log(x),3)+7.266e-03*pow(log(x),4)"},}},

  {"trkEffNtrk1m3",          {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk1m3"},      {"appliesTo","Resp"}, {"form", "+6.020e+00-5.647e+00*log(x)+1.649e+00*pow(log(x),2)-1.955e-01*pow(log(x),3)+8.229e-03*pow(log(x),4)"},}},
  {"trkEffNtrk1m3_chf",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk1m3"},      {"appliesTo","chf"},  {"form", "+4.784e+00-5.104e+00*log(x)+1.600e+00*pow(log(x),2)-1.994e-01*pow(log(x),3)+8.710e-03*pow(log(x),4)"},}},
  {"trkEffNtrk1m3_nhf",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk1m3"},      {"appliesTo","nhf"},  {"form", "-2.550e+00+2.776e+00*log(x)-8.980e-01*pow(log(x),2)+1.156e-01*pow(log(x),3)-5.216e-03*pow(log(x),4)"},}},
  {"trkEffNtrk1m3_nef",      {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk1m3"},      {"appliesTo","nef"},  {"form", "-3.124e+00+3.086e+00*log(x)-9.373e-01*pow(log(x),2)+1.149e-01*pow(log(x),3)-4.982e-03*pow(log(x),4)"},}},

  {"trkEffNtrk2ToInfm3",     {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk2ToInfm3"}, {"appliesTo","Resp"}, {"form", "+3.457e+00-1.972e+00*log(x)+2.966e-01*pow(log(x),2)-9.399e-03*pow(log(x),3)-4.569e-04*pow(log(x),4)"},}},
  {"trkEffNtrk2ToInfm3_chf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk2ToInfm3"}, {"appliesTo","chf"},  {"form", "-8.492e+00+7.806e+00*log(x)-2.447e+00*pow(log(x),2)+2.963e-01*pow(log(x),3)-1.222e-02*pow(log(x),4)"},}},
  {"trkEffNtrk2ToInfm3_nhf", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk2ToInfm3"}, {"appliesTo","nhf"},  {"form", "-4.057e+00+2.813e+00*log(x)-6.923e-01*pow(log(x),2)+7.966e-02*pow(log(x),3)-3.445e-03*pow(log(x),4)"},}},
  {"trkEffNtrk2ToInfm3_nef", {{"ispositive", "0"}, {"initial", "0"}, {"freeze", "0"}, {"type","trkEffNtrk2ToInfm3"}, {"appliesTo","nef"},  {"form", "+1.414e+01-1.197e+01*log(x)+3.551e+00*pow(log(x),2)-4.294e-01*pow(log(x),3)+1.815e-02*pow(log(x),4)"},}},

};
// clang-format on