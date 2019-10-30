void GetAllInfo(){

  // Survey Function to Get Run list from SQL database
  // Author : Tao Ye

  TSQLResult* res;
  TSQLRow *row;
  cout << " -- Getting Slug Number from RCDB -- " << endl;
  cout << " -- Connecting to RCDB -- " << endl;
  TSQLServer *rcdb = TSQLServer::Connect("mysql://hallcdb.jlab.org:3306/a-rcdb","rcdb","");
  cout << " -- ServerInfo: " << rcdb->ServerInfo() << endl;
  cout << " -- Host : " << rcdb->GetHost() << endl;
  cout << " -- Query DataBase " << endl;

  TString query[]={"SELECT t1.run_number,t3.int_value,t4.text_value, t1.text_value , t2.text_value,CASE WHEN tflag.text_value is NULL THEN 'NoFlag' ELSE tflag.text_value END AS flag_res ,t6.int_value, CASE WHEN tcomt.text_value is NULL THEN '- No comment -' ELSE tcomt.text_value END",
  		   " FROM `a-rcdb`.conditions as t1 INNER JOIN `a-rcdb`.condition_types c1 on c1.id= t1.condition_type_id AND c1.name='ihwp' ",
  		   ", `a-rcdb`.conditions as t2 INNER JOIN `a-rcdb`.condition_types c2 on c2.id=t2.condition_type_id AND c2.name='flip_state' ",
  		   ", `a-rcdb`.conditions as t3 INNER JOIN `a-rcdb`.condition_types c3 on c3.id=t3.condition_type_id AND c3.name='slug' AND t3.int_value BETWEEN 0 AND 94 ",
  		   ", `a-rcdb`.conditions as t4 INNER JOIN `a-rcdb`.condition_types c4 on c4.id=t4.condition_type_id AND c4.name='run_type' ",
  		   " LEFT JOIN ",
		   "( SELECT run_number, text_value ",
		   " FROM `a-rcdb`.conditions t5 INNER JOIN `a-rcdb`.condition_types c5 on c5.id=t5.condition_type_id AND c5.name='run_flag' ) AS tflag ",
		   " ON tflag.run_number = t4.run_number ",
  		   ", `a-rcdb`.conditions as t6 INNER JOIN `a-rcdb`.condition_types c6 on c6.id=t6.condition_type_id AND c6.name='arm_flag' ",
  		   " LEFT JOIN ",
		   "( SELECT run_number, text_value ",
		   " FROM `a-rcdb`.conditions t7 INNER JOIN `a-rcdb`.condition_types c7 on c7.id=t7.condition_type_id AND c7.name='wac_comment' ) AS tcomt ",
		   " ON tcomt.run_number = t6.run_number ",
  		   " WHERE t1.run_number=t3.run_number AND t1.run_number=t2.run_number AND t1.run_number=t4.run_number AND t1.run_number=t6.run_number AND t1.run_number<4981",
		   " ORDER BY t1.run_number ASC"};

  TString listName = Form("./survey.list");
  
  int n=sizeof(query)/sizeof(*query);
  TString cmd;
  for(int  i=0;i<n;i++){
    cmd+=query[i];
  }
  res = rcdb->Query(cmd);

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

  FILE *runlist = fopen(listName.Data(),"w");
  for(int irow=0;irow<nRows;irow++){
    row=res->Next();
    for(int j=0; j< nFields; j++){
      cout << "\t" << row->GetField(j) ;
      TString row_ret = TString(row->GetField(j));
      row_ret.ReplaceAll("\n",";");
      fprintf(runlist,"%s",row_ret.Data());
      if(j!=nFields-1)
	fprintf(runlist,",");
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
