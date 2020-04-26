TString generate_cut(vector<Double_t> fEventCounterArray){
  TString fret;
  Int_t length = fEventCounterArray.size();
  if(length==0)
    return "1==0";
  
  Int_t start_pt= fEventCounterArray[0];
  Int_t end_pt = -1;
  vector<pair<Int_t,Int_t> > fRangeArray;
  for(int i=1;i<length;i++){
    if(fEventCounterArray[i]-fEventCounterArray[i-1]>1){
      end_pt = fEventCounterArray[i-1];
      fRangeArray.push_back(make_pair(start_pt,end_pt));
      start_pt = fEventCounterArray[i];
    }

    if(i==length-1){
      end_pt = fEventCounterArray[i];
      fRangeArray.push_back(make_pair(start_pt,end_pt));
    }
  }
  auto iter = fRangeArray.begin();
  while(iter!=fRangeArray.end()){
    if( (*iter).first == (*iter).second){
      iter++;
      continue;
    }
    if(fret=="")
      fret += Form("(pattern_number>=%d && pattern_number<=%d )",
		   (*iter).first, (*iter).second);
    else
      fret += Form("|| (pattern_number>=%d && pattern_number<=%d )",
		   (*iter).first, (*iter).second);
    
    iter++;
  }
  
  return fret;
}
