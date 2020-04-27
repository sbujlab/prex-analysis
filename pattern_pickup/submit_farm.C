void submit_farm(){

  // FILE *output = fopen("by_block_all_slug.xml","w");
  FILE *output = fopen("by_block_all_slug.xml","w");
  fprintf(output,"<Request>\n");
  fprintf(output,"  <Email email=\"tao@jlab.org\" request=\"false\" job=\"true\" />\n");
  fprintf(output,"  <Project name=\"prex\"/> \n");
  fprintf(output,"  <Track name=\"one_pass\"/> \n");
  fprintf(output,"  <Name name=\"by_block_all_slug\"/> \n");
  fprintf(output,"  <OS name=\"centos7\"/> \n");
  fprintf(output,"  <Memory space=\"1000\" unit=\"MB\" /> \n");

  for(int islug=1;islug<=94;islug++){
    fprintf(output,"  <Job>\n");
    fprintf(output,"    <Command><![CDATA[\n");
    fprintf(output,"      setenv QW_ROOTFILES /lustre/expphy/volatile/halla/parity/japanOutput/ \n");
    fprintf(output,"      cd /u/home/tao/work/pattern-pickup/ \n");
    // fprintf(output,"      root -b -q 'GetAsymmetryByBlockFromCut_quick2.C(%d)'  \n",islug);
    fprintf(output,"      root -b -q 'GetAsymmetryByBlock_new.C(%d)'  \n",islug);
    fprintf(output,"    ]]></Command>\n");
    fprintf(output,"  <Stdout dest=\"/u/home/tao/work/pattern-pickup/farm_log/farm_log_slug%d_byBlock_new.out\"/> \n ",islug );
    fprintf(output,"  <Stderr dest=\"/u/home/tao/work/pattern-pickup/farm_log/farm_log_slug%d_byBlock_new.err\"/> \n",islug );
    fprintf(output,"  </Job>\n");
    fprintf(output,"\n");
  }

  fprintf(output,"</Request>\n");
  fclose(output);
}
