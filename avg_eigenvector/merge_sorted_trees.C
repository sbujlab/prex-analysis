void merge_sorted_trees(TString type = "", TString suffix = "") { // Type holds "part_avg" for 2nd round ana outputs
  //TFile *_file0 = TFile::Open("rootfiles/mini_eigen_reg_allbpms_sorted.root");
  //TFile *_file1 = TFile::Open("rootfiles/mini_eigen_reg_5bpms_sorted.root");
  //TFile *_file2 = TFile::Open(Form("dataRootfiles/rcdb_eigenvectors.root");
  //
  //FIXME Default analysis outputs, Feb 24th
  if (type == "") {
    TFile *ofile = TFile::Open(Form("rootfiles/rcdb_eigenvectors_sorted%s.root",suffix.Data()),"recreate");
    TChain* mini_eigen_reg_allbpms_sorted = new TChain(Form("mini_eigen_reg_allbpms_sorted"));
    mini_eigen_reg_allbpms_sorted->Add(Form("rootfiles/mini_eigen_reg_allbpms_sorted%s.root",suffix.Data()));
    TChain* mini_eigen_reg_5bpms_sorted = new TChain(Form("mini_eigen_reg_5bpms_sorted"));
    mini_eigen_reg_5bpms_sorted->Add(Form("rootfiles/mini_eigen_reg_5bpms_sorted%s.root",suffix.Data()));
    TChain* mini = new TChain("mini");
    /*mini->Add(Form("rootfiles/rcdb_eigenvectors%s.root",suffix.Data()));
    TChain* mini_eigen_reg_allbpms = new TChain("mini_eigen_reg_allbpms");
    mini_eigen_reg_allbpms->Add(Form("rootfiles/rcdb_eigenvectors%s.root",suffix.Data()));
    TChain* mini_eigen_reg_allbpms_tr = new TChain("mini_eigen_reg_allbpms_tr");
    mini_eigen_reg_allbpms_tr->Add(Form("rootfiles/rcdb_eigenvectors%s.root",suffix.Data()));
    TChain* mini_eigen_reg_5bpms = new TChain("mini_eigen_reg_5bpms");
    mini_eigen_reg_5bpms->Add(Form("rootfiles/rcdb_eigenvectors%s.root",suffix.Data()));
    */
    mini->Add(Form("rootfiles/rcdb_eigenvectors.root"));
    TChain* mini_eigen_reg_allbpms = new TChain("mini_eigen_reg_allbpms");
    mini_eigen_reg_allbpms->Add(Form("rootfiles/rcdb_eigenvectors.root"));
    TChain* mini_eigen_reg_allbpms_tr = new TChain("mini_eigen_reg_allbpms_tr");
    mini_eigen_reg_allbpms_tr->Add(Form("rootfiles/rcdb_eigenvectors.root"));
    TChain* mini_eigen_reg_5bpms = new TChain("mini_eigen_reg_5bpms");
    mini_eigen_reg_5bpms->Add(Form("rootfiles/rcdb_eigenvectors.root"));
    ofile->cd();
    mini_eigen_reg_allbpms_tr->CloneTree()->Write("",TObject::kOverwrite);
    mini_eigen_reg_allbpms->CloneTree()->Write("",TObject::kOverwrite);
    mini_eigen_reg_allbpms_sorted->CloneTree()->Write("",TObject::kOverwrite);
    mini_eigen_reg_5bpms_sorted->CloneTree()->Write("",TObject::kOverwrite);
    mini_eigen_reg_5bpms->CloneTree()->Write("",TObject::kOverwrite);
    mini->CloneTree()->Write("",TObject::kOverwrite);
    ofile->Close();
  }
  else {

    // FIXME March 16th part-wise eigen vector outputs
    TFile *ofile = TFile::Open(Form("rootfiles/rcdb_eigenvectors_sorted%s%s.root",type.Data(),suffix.Data()),"recreate");
    //TFile *ofile = TFile::Open(Form("rootfiles/rcdb_eigenvectors_sorted_part_avg%s.root",suffix.Data()),"recreate");

    TChain* mini = new TChain("mini");
    mini->Add(Form("rootfiles/rcdb_eigenvectors%s.root",suffix.Data()));

    TChain* mini_eigen_reg_5bpms_part_avg = new TChain(Form("mini_eigen_reg_5bpms_part_avg"));
    mini_eigen_reg_5bpms_part_avg->Add(Form("rootfiles/rcdb_eigenvectors%s.root",suffix.Data()));

    TChain* mini_reference_eigen_reg_5bpms = new TChain(Form("mini_reference_eigen_reg_5bpms%s",suffix.Data()));
    mini_reference_eigen_reg_5bpms->Add(Form("rootfiles/rcdb_eigenvectors%s.root",suffix.Data()));
    TChain* mini_reference_eigen_reg_5bpms_sorted = new TChain(Form("mini_reference_eigen_reg_5bpms_sorted%s",suffix.Data()));
    mini_reference_eigen_reg_5bpms_sorted->Add(Form("rootfiles/mini_reference_eigen_reg_5bpms_sorted%s.root",suffix.Data()));

    ofile->cd();

    mini->CloneTree()->Write("",TObject::kOverwrite);
    mini_eigen_reg_5bpms_part_avg->CloneTree()->Write("",TObject::kOverwrite);
    mini_reference_eigen_reg_5bpms->CloneTree()->Write("",TObject::kOverwrite);
    mini_reference_eigen_reg_5bpms_sorted->CloneTree()->Write("",TObject::kOverwrite);

    ofile->Close();
  }
}
