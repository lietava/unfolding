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
        TF1 *fun = new TF1("fun","(sin(x))",0,2*TMath::Pi());
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
int test11(){
	//gSystem->Load("/home/rl/RooUnfold/libRooUnfol");	
	Int_t Ndata = 10000; 
	Double_t* datagen = new Double_t[Ndata];
	Double_t* datarec = new Double_t[Ndata];
        Double_t* datagen_sin = new Double_t[Ndata];
        Double_t* datarec_sin = new Double_t[Ndata];
	//
        generate_data(Ndata,datagen);
        generate_data_sin(Ndata,datagen_sin);
	reconstruct_data(Ndata,datagen,datarec);
        reconstruct_data(Ndata,datagen_sin,datarec_sin);
	// Create bins
        //std::vector<Double_t> vbins = {0,1.14,2,3,4,5,2*TMath::Pi()};
        std::vector<Double_t> vbins = uniformbins(11);
	//
	Int_t nbins = vbins.size()-1;
	Double_t *bins = new Double_t[nbins+1];
	for(int i = 0; i < (nbins+1); i++) bins[i] = vbins[i];
        std:cout << "N data bins:" << nbins << std::endl;
	// Create histos
        TH1D *hdatagen = new TH1D("data gen","data generated",nbins,bins);
        TH1D *hdatagen_rec = new TH1D("data gen rec","data generated&reconstructed",nbins,bins);
        TH1D *hdatarec = new TH1D("data rec","data reconstructed",nbins,bins);
        TH1D *hres = new TH1D("resolution","resolution",20,-2,2);
        //Double_t* datagenh = datagen_sin;
        //Double_t* datarech = datarec_sin;
        Double_t* datagenh = datagen;
        Double_t* datarech = datarec;
	for(int i = 0; i < Ndata; i++) {
                hdatagen->Fill(datagenh[i]);
                if(datarech[i] >= 0) {
                        hdatarec->Fill(datarech[i]);
                        hdatagen_rec-> Fill(datagenh[i]);
                        hres->Fill(datagenh[i]-datarech[i]);
		}
	}
	// Response
	RooUnfoldResponse response(hdatarec,hdatagen);
        // choose data set
        Double_t* datagenr = datagen;
        Double_t* datarecr = datarec;
	for(int i  = 0; i < Ndata; i++) {
		//std::cout << i << " response" << " gen:" << datagen[i] << " rec:" << datarec[i] << std::endl;
                if(datarecr[i] >= 0)
                  response.Fill(datarecr[i],datagenr[i]);
		else
                  response.Miss(datagenr[i]);
	}
        // Generated folded
	TH1D* hfolded = (TH1D*)response.ApplyToTruth(hdatagen);
	//
        // Reconstructed Unfolded
        RooUnfoldBayes unfold(&response, hdatarec, 10);
        //RooUnfoldBayes unfold("bayes","bayes");
        unfold.SetVerbose(1);
        //unfold.SetResponse(&response);
        //unfold.SetMeasured(hdatarec);
        //unfold.SetIterations(10);
        //unfold.SetSmoothing(0);
        unfold.SetPriors(hdatagen);
        //unfold.Unfold();
        TH1D* hunfold = (TH1D*) unfold.Hreco();
        //unfold.SetPriors(hdatarec);

        //return 0;
        Double_t chi2 = calcChi2(hdatagen,hunfold,nbins);
        std::cout << "chi2 gen-unfold:" << chi2 << std::endl;
	// Folding back
	TH1D* hrefold = (TH1D*)response.ApplyToTruth(hunfold);
        chi2 = calcChi2(hdatarec,hrefold,nbins);
        std::cout << "chi2 rec-refold:" << chi2 << std::endl;

	//
	TFile *myfile = new TFile("unfold.root", "RECREATE");
	hdatagen->Write();
	hdatarec->Write();
        hunfold->Write("unfolded: R-1*meas");
        hfolded->Write("folded: R*gen");
        hrefold->Write("hrefold: R*unfolded");
        hdatagen_rec->Write("Data rec");
        hres->Write();
	response.Hresponse()->Write("Response");
	myfile->Close();
	return 0;
}
