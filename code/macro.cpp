{
gROOT->ProcessLine("gSystem->Load(\"$LIB/libTKinFitter.so\")");
gROOT->ProcessLine(".L events.cpp+");
gROOT->ProcessLine("events()");
}
