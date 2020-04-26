void GetRunList(Int_t run_number=4040){

  // Experimenting Function to Get slug number based on run number 
  // Author : Tao Ye

  TSQLResult* res;
  TSQLRow *row;
  cout << " -- Getting Slug Number from RCDB -- " << endl;
  cout << " -- Connecting to RCDB -- " << endl;
  TSQLServer *rcdb = TSQLServer::Connect("mysql://hallcdb.jlab.org:3306/a-rcdb","rcdb","");
  cout << " -- ServerInfo: " << rcdb->ServerInfo() << endl;
  cout << " -- Host : " << rcdb->GetHost() << endl;
  rcdb->GetTablesList()->Print();
  rcdb->GetTableInfo("runs")->GetColumns()->Print();
  rcdb->GetTableInfo("conditions")->GetColumns()->Print();
  rcdb->GetTableInfo("condition_types")->GetColumns()->Print();
  cout << " -- Query DataBase " << endl;
  TString select_q ="SELECT run_number,int_value "; 
  TString from_q =  "FROM `a-rcdb`.conditions INNER JOIN `a-rcdb`.condition_types on  condition_types.id=condition_type_id ";
  TString order_q = "ORDER BY run_number ASC ";
  TString where_q =  "WHERE run_number in (SELECT run_number FROM `a-rcdb`.conditions INNER JOIN `a-rcdb`.condition_types on condition_types.id=condition_type_id and text_value='Production' and name='run_type' ) and int_value between 0 and 94 and name='slug' ";

  res = rcdb->Query(select_q + from_q + where_q + order_q);

  if(res==NULL){
    cout << " -- Failed to Query " << endl;
    cout << " -- Bye-bye! " << endl;
    delete res;
    cout << " -- Closing Connection to RCDB " << endl;
    rcdb->Close();
    delete rcdb;
    return -1;
  }

  int nFields =res->GetFieldCount();
  int nRows = res->GetRowCount();
  if(nRows==0){
    cout << " -- Failed to load Results " << endl;
    cout << " -- Bye-bye! " << endl;
    delete res;
    cout << " -- Closing Connection " << endl;
    rcdb->Close();
    delete rcdb;
    return ;
  }
  cout << " ----------------------------- " << endl;
  FILE *runlist = fopen("all_production.list","w");
  for(int irow=0;irow<nRows;irow++){
    row=res->Next();
    for(int j=0; j< nFields; j++){
      cout << "\t" << row->GetField(j) ;
      fprintf(runlist,"%s\t",row->GetField(j));
    }
    fprintf(runlist,"\n");
    cout << endl;
  }
  delete row;
  cout << " ----------------------------- " << endl;
  cout << " -- Closing Connection " << endl;
  rcdb->Close();
  delete res;
  delete rcdb;
}
