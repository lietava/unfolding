Double_t sigma = 2;
TRandom3 rnd;
//std::vector<Double_t> bins = {0,1,3.14,5,2*TMath::Pi()}; 
int generate_data(int Ndata,Double_t* gen) {
	rnd.RndmArray(Ndata,gen);
	for(int i =0; i<Ndata; i++) gen[i] = gen[i]*2*TMath::Pi();
	return 0;
}
//
int generate_data_sin(int Ndata,Double_t* gen) {
        TF1 *fun = new TF1("fun","sin(x)",0,2*TMath::Pi());
        for(int i =0; i<Ndata; i++) gen[i] = fun->GetRandom();
        return 0;
}
// parabolic efficiency
Double_t efficiency(Double_t phi) {
	Double_t A=2*TMath::Pi();
	return 4*phi*(A-phi)/A/A;
}
int reconstruct_data(int Ndata,Double_t* datagen, Double_t* datarec) {	
	Double_t* rr = new Double_t[Ndata];
	rnd.RndmArray(Ndata,rr);
	for(int i =0; i < Ndata; i++) {
		Double_t eff = efficiency(datagen[i]);
		if(eff < rr[i]) {
			datarec[i] = -1;
		} else {
			Double_t phimeas = rnd.Gaus(datagen[i], sigma);
			//if(phimeas < 0) phimeas += 2*TMath::Pi();
			//if(phimeas < 0) phimeas = rnd.Rndm()*2*TMath::Pi();
			//if(phimeas > 2*TMath::Pi()) phimeas = rnd.Rndm()*2*TMath::Pi();
			if((phimeas < 0) || (phimeas > 2*TMath::Pi())) phimeas = remainder(phimeas,2*TMath::Pi());
			if(phimeas < 0) phimeas += 2*TMath::Pi();
			datarec[i] = phimeas;
		}
                //if(datagen[i] > 5.) std::cout << i << " eff:" << eff << " "<< rr[i] << " gen:" << datagen[i] << " rec:" << datarec[i] << std::endl;
	}
	return 0;
}
std::vector<Double_t> uniformbins(int nbins) {
	std::vector<Double_t> vbins;
	int nb = nbins;
	double del = 2*TMath::Pi()/(nb-1);
	for(int i =0; i < nb; i++) vbins.push_back(i*del);
	return vbins;
}
Double_t calcChi2(TH1* h1, TH1* h2, int nbins) {
    Int_t nb = h1->GetNbinsX();
    std::cout << "chi2 Nbins:" << nb << std::endl;
    Double_t chi2 = 0;
    for(int i = 1; i < nbins+1; i++) {
        Double_t b1 =h1->GetBinContent(i);
        Double_t b2 =h2->GetBinContent(i);
        //std::cout << "b1:" << b1 << " b2:" << b2 << std::endl;
        chi2 += (b1-b2) *(b1-b2)/b1;
    }
    return chi2;
}
RooUnfoldResponse* getResponseFromFile() {
    TFile *file = new TFile("download/UnfoldedClosureHe_10.root");
    RooUnfoldResponse* response = (RooUnfoldResponse*) file->Get("response;1");
    return response;
}
int test2(){
	//gSystem->Load("/home/rl/RooUnfold/libRooUnfol");	
        //
        RooUnfoldResponse * res = getResponseFromFile();
        //res->Hresponse()->Draw();
        //return 0;
	Int_t Ndata = 10000; 
	Double_t* datagen = new Double_t[Ndata];
        //Double_t* datarec = new Double_t[Ndata];
        Double_t* datagen2 = new Double_t[Ndata];
        //Double_t* datarec2 = new Double_t[Ndata];
	//
        generate_data(Ndata,datagen);
        generate_data_sin(Ndata,datagen2);

        //
        
	// Create histos
        TH1D *hdatagen = (TH1D*)res->Htruth()->Clone();
        Int_t nbins = hdatagen->GetNbinsX();
        hdatagen->Reset("ICESM");
        //TH1D *hdatagen_rec = new TH1D("data gen rec","data generated&reconstructed",nbins,bins);
        //TH1D *hdatarec = res->Hmeasured()->Clone();
        //TH1D* hfolded = (TH1D*)response.ApplyToTruth(hdatagen);
        Double_t* datagenh = datagen2;
        //Double_t* datarech = datarec;
        for(int i = 0; i < Ndata; i++) {
                hdatagen->Fill(datagenh[i]);
        }
        // Response
        RooUnfoldResponse* response = res;;
        // Reconstruct
        TH1D* hdatarec = (TH1D*)response->ApplyToTruth(hdatagen);

        //
        // Reconstructed Unfolded
        RooUnfoldBayes unfold(response, hdatarec, 10);
	TH1D* hunfold = (TH1D*) unfold.Hreco();
        Double_t chi2 = calcChi2(hdatagen,hunfold,nbins);
        std::cout << "chi2 gen-unfold:" << chi2 << std::endl;
	// Folding back
        TH1D* hrefold = (TH1D*)response->ApplyToTruth(hunfold);
        chi2 = calcChi2(hdatarec,hrefold,nbins);
        std::cout << "chi2 rec-refold:" << chi2 << std::endl;

	//
	TFile *myfile = new TFile("unfold.root", "RECREATE");
        hdatagen->Write("generated");
	hdatarec->Write();
	hunfold->Write("unfolded");
        //hfolded->Write("folded");
	hrefold->Write("hrefold");
        //hdatagen_rec->Write();
        //hres->Write();
        response->Hresponse()->Write("Response");
	myfile->Close();
	return 0;
}
