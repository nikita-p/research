void macro2(string file, string key){
    gROOT->ProcessLine("gSystem->Load(\"$LIB/libTKinFitter.so\")");
    gROOT->ProcessLine(".L events.cpp+");
    string s = "events_cores(\""+file+"\", \"" + key + "\")";
//     cout << s << endl;
    gROOT->ProcessLine(s.c_str());
    gROOT->ProcessLine(".q");
    return;
}
