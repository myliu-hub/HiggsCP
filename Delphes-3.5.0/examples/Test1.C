void Test1()
{

   TChain chain("Delphes");
   ifstream fp;
   fp.open("/afs/ihep.ac.cn/users/m/myliu/scratchfs/Delphes-3.5.0/all_txt/Pcc.txt");
   char buf[200];
   int a =0;
   while(fp.getline(buf,200)){
cout << "a= " << a << endl;
     chain.Add(buf);
     a++;
   }
   TFile * newfile =TFile::Open("/afs/ihep.ac.cn/users/m/myliu/scratchfs/Delphes-3.5.0/Data/result/test_Pcc.root","recreate");
   TTree * tr1 =chain.CloneTree();
   newfile->cd();
   tr1->Write();
   newfile->Write();
   newfile->Close();
}
