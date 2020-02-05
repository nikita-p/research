#include <iostream>
#include <TComplex.h>
#include <TF1.h>

using namespace std;

double fas(double TwoE)
{
    const int K_max = 4;

    bool Calc = true;

    Double_t A[K_max] = {-6.12394622, 25.0341405, -34.1311022, 15.5413717};
    Double_t B[K_max] = {5.29354148, -7.90990714, -2.26007613, 5.21453902};
    Double_t C[K_max] = {-0.56436115415016, 2.69953793, -4.32966739, 2.33116866};
    Double_t D[K_max] = {-0.0548334238, 0.31600391, -0.609523718, 0.393667808};

    Double_t A1[K_max] = {-4.91624401, 19.8655606, -26.9136128, 12.2412286};
    Double_t D1[K_max] = {-0.00794774189, 0.0522269164, -0.114526409, 0.0838126536};

    Double_t LM[6] = {1.1, 0.875, 0.75, 0.62, 0.52, 0.46};

    Double_t POL, TEMP1, TEMP2;

    TEMP2 = D[1] * LM[4] + D[2] * LM[4] * LM[4] + D[3] * LM[4] * LM[4] * LM[4];
    TEMP1 = D1[0] + D1[1] * LM[4] + D1[2] * LM[4] * LM[4] + D1[3] * LM[4] * LM[4] * LM[4];
    D[0] = TEMP1 - TEMP2;

    TEMP2 = C[1] * LM[3] + C[2] * LM[3] * LM[3] + C[3] * LM[3] * LM[3] * LM[3];
    TEMP1 = D[0] + D[1] * LM[3] + D[2] * LM[3] * LM[3] + D[3] * LM[3] * LM[3] * LM[3];
    C[0] = TEMP1 - TEMP2;

    TEMP2 = A1[1] * LM[2] + A1[2] * LM[2] * LM[2] + A1[3] * LM[2] * LM[2] * LM[2];
    TEMP1 = C[0] + C[1] * LM[2] + C[2] * LM[2] * LM[2] + C[3] * LM[2] * LM[2] * LM[2];
    A1[0] = TEMP1 - TEMP2;

    TEMP1 = A1[0] + A1[1] * LM[1] + A1[2] * LM[1] * LM[1] + A1[3] * LM[1] * LM[1] * LM[1];
    TEMP2 = A[1] * LM[1] + A[2] * LM[1] * LM[1] + A[3] * LM[1] * LM[1] * LM[1];
    A[0] = TEMP1 - TEMP2;

    TEMP1 = A[0] + A[1] * LM[0] + A[2] * LM[0] * LM[0] + A[3] * LM[0] * LM[0] * LM[0];
    TEMP2 = B[1] * LM[0] + B[2] * LM[0] * LM[0] + B[3] * LM[0] * LM[0] * LM[0];
    B[0] = TEMP1 - TEMP2;

    Calc = false;

    POL = 0.0;
    if (TwoE >= LM[0])
        for (Int_t i = 0; i < K_max; i++)
            POL += TMath::Power(TwoE, (Double_t)i) * B[i];

    if (TwoE >= LM[1] && TwoE < LM[0])
        for (Int_t i = 0; i < K_max; i++)
            POL += TMath::Power(TwoE, (Double_t)i) * A[i];

    if (TwoE >= LM[2] && TwoE < LM[1])
        for (Int_t i = 0; i < K_max; i++)
            POL += TMath::Power(TwoE, (Double_t)i) * A1[i];

    if (TwoE >= LM[3] && TwoE < LM[2])
        for (Int_t i = 0; i < K_max; i++)
            POL += TMath::Power(TwoE, (Double_t)i) * C[i];

    if (TwoE >= LM[4] && TwoE < LM[3])
        for (Int_t i = 0; i < K_max; i++)
            POL += TMath::Power(TwoE, (Double_t)i) * D[i];

    if (TwoE >= LM[5] && TwoE < LM[4])
        for (Int_t i = 0; i < K_max; i++)
            POL += TMath::Power(TwoE, (Double_t)i) * D1[i];
    if (TwoE < LM[5])
        POL = D1[0] + D1[1] * LM[5] + D1[2] * LM[5] * LM[5] + D1[3] * LM[5] * LM[5] * LM[5];

    double Fval = (POL / 0.393728 * (0.00749 / 0.0361478));

    return Fval;
}

class MDVM
{

    static const double ALPHA;
    static const double C;

    static const double m_phi;
    static const double m_rho;
    static const double m_omg;

    static const double m_k0;
    static const double m_pi;
    static const double m_kc;
    static const double m_pi0;

    static const double w0_phi;
    static const double w0_rho;
    static const double w0_omg;
public:
    static double BETA(double s, double M_K);
    //бета
    static double PV2(double s, double M, double Mn);
    //PV2: фазовый объём (почти: он типа нормирован) распада на 2 одинаковые частицы, s: энергия**2, M: масса векторного мезона, W0: ширина распада,  Mn: масса частицы в распаде
    static double PV3(double s, double M, double m1, double m2, double m3);
    //PV3: фазовый объём (почти: он типа нормирован) распада на 3 частицы
    static double PVG(double s, double M, double Mn);
    //PVG: распад частицы M на фотон и частицу Mn

    static double WOmg(double s, double W0, double MX);

    static double WRhoX(double s, double W0, double MX)
    {
        return W0 * PV2(s, MX, m_pi);
    }
    static double WOmgX(double s, double W0, double MX)
    {
        return W0 * PV3(s, MX, m_pi, m_pi, m_pi0);
    }
    static double WPhiX(double s, double W0, double MX)
    {
        return W0 * PV2(s, MX, m_kc);
    }

    static TComplex BW(double, double, double, double (*WX)(double, double, double));
    //функция Брейта-Вигнера
    static TComplex BW_Rho(double s)
    {
        return BW(s, m_rho, w0_rho, WRhoX);
    }
    static TComplex BW_Omg(double s)
    {
        return BW(s, m_omg, w0_omg, WOmg);
    }
    static TComplex BW_Phi(double s)
    {
        return BW(s, m_phi, w0_phi, WPhi);
    }

    static TComplex BW_Rho1(double s)
    { // 1465, 25; 400, 60;
        return BW(s, 1490, 340, WRhoX);
    }
    static TComplex BW_Omg1(double s)
    {
        return BW(s, 1420, 220, WOmgX);
    }
    static TComplex BW_Rho2(double s)
    {
        return BW(s, 1574, 234, WRhoX);
    }
    static TComplex BW_Omg2(double s)
    {
        return BW(s, 1688, 350, WOmgX);
    }
    static TComplex BW_Phi1(double s)
    {
        return BW(s, 1673, 182, WPhiX);
    }
    static TComplex BW_Rho3(double s)
    { //unused
        return BW(s, 1720, 250, WRhoX);
    }
    static TComplex BW_Rho4(double s)
    { //unused
        return BW(s, 1880, 160, WRhoX);
    }
    static TComplex BW_Rho5(double s)
    {
        return BW(s, 2134, 343, WRhoX);
    }
    static TComplex BW_Phi2(double s)
    {
        return BW(s, 2198, 71, WPhiX);
    }

//public:
    static double WPhi(double s, double W0, double MX);
    static TComplex F0(double *x, double *par, bool mode);
    //формфактор, нулевое приближение
    //mode: 0 - short/long; 1 - charged;
    static TComplex F1(double *x, double *par, bool mode);
    //формфактор с учётом omega(1400)
    static double Cross_Section(double *x, double *par, bool mode);
    static double Cross_Section_Neutral(double *x, double *par)
    {
        return Cross_Section(x, par, 0);
    }
    static double Cross_Section_Charged(double *x, double *par)
    {
        return Cross_Section(x, par, 1);
    }

    static TF1 *Cross_Section(bool mode);
};

/*___________________________________________________________*/

const double MDVM::C = 0.389379292E12; //(MeV)^2 * nb
const double MDVM::ALPHA = 7.297E-3;

const double MDVM::m_rho = 775.26;
const double MDVM::m_omg = 782.65;
const double MDVM::m_phi = 1019.464; //1;

const double MDVM::m_k0 = 497.6;
const double MDVM::m_kc = 493.677;
const double MDVM::m_pi = 139.57;
const double MDVM::m_pi0 = 135.;

const double MDVM::w0_phi = 4.247; //9;
const double MDVM::w0_rho = 147.8;
const double MDVM::w0_omg = 8.49;

/*__________________________________________________________*/

double MDVM::BETA(double s, double M_K)
{ //бета
    double E = TMath::Sqrt(s) / 2.;
    double P = TMath::Sqrt(E * E - M_K * M_K);
    return P / E;
}

double MDVM::PV2(double s, double M, double Mn)
{
    if (s / 4. < Mn * Mn)
        return 0;
    double E = sqrt(s) / 2;
    double w = (M * M / s) * pow((s / 4. - Mn * Mn) / ((M * M) / 4. - Mn * Mn), 3 / 2.);

    return w;
}

/*double MDVM::PV3(double s, double MX, double m1, double m2, double m3){
        double pv = (pow(TMath::Pi(),3)/2.)*( pow(m1*m2*m3,1/2.)*pow(sqrt(s) - m1 - m2 - m3,2)/pow(m1 + m2 + m3,3/2.) );
        double pv0 = (pow(TMath::Pi(),3)/2.)*( pow(m1*m2*m3,1/2.)*pow(sqrt(MX*MX) - m1 - m2 - m3,2)/pow(m1 + m2 + m3,3/2.) );
        return pv/pv0;
    }*/

double MDVM::PV3(double s, double MX, double m1, double m2, double m3)
{
    double pv = fas(sqrt(s) / 1000);
    return pv / fas(MX / 1000);
}

double MDVM::PVG(double s, double MX, double Mn)
{
    double pv = pow((s - Mn * Mn) / (2 * sqrt(s)), 3);
    double pv0 = pow((MX * MX - Mn * Mn) / (2 * sqrt(MX * MX)), 3);
    return pv / pv0;
}

double MDVM::WPhi(double s, double W0, double MX)
{
    double Br_KC = 0.492;
    double Br_KN = 0.34;
    double Br_3Pi = 0.1524;
    double Br_EG = 0.01303;
    double ost = 1 - Br_KC - Br_KN - Br_3Pi - Br_EG;

    double m_eta = 547.862;

    double W = W0 * ((Br_KC + ost) * PV2(s, MX, m_kc) + Br_KN * PV2(s, MX, m_k0) + Br_3Pi * PV3(s, MX, m_pi, m_pi, m_pi0) + Br_EG * PVG(s, MX, m_eta));
    return W;
}

double MDVM::WOmg(double s, double W0, double MX)
{
    double Br_3Pi = 0.892;
    double Br_Pi0G = 0.084;
    double Br_2Pi = 0.0153;
    double ost = 1 - Br_Pi0G - Br_2Pi;

    double W = W0 * ((Br_3Pi + ost) * PV3(s, MX, m_pi, m_pi, m_pi0) + Br_Pi0G * PVG(s, MX, m_pi0) + Br_2Pi * PV2(s, MX, m_pi));
    return W;
}

TComplex MDVM::BW(double s, double MX, double WX0, double (*WX)(double, double, double))
{
    TComplex I(0, 1);
    TComplex bw = pow(MX, 2) / (pow(MX, 2) - s - I * sqrt(s) * WX(s, WX0, MX));
    return bw;
}

TComplex MDVM::F0(double *x, double *par, bool mode)
{
    double n = par[10]; //1.026;//1.027;
    double s = TMath::Power(x[0] * 1E3, 2);
    double CR = par[0];
    double CO = par[1];
    double CP = par[2];

    double KR = mode ? CR / 2. : -CR / 2.;
    double KO = CO / 6.;
    double KP = mode ? CP / 3. : n * CP / 3.;

    TComplex F = KR * BW_Rho(s) + KO * BW_Omg(s) + KP * BW_Phi(s);
    return F;
}

TComplex MDVM::F1(double *x, double *par, bool mode)
{
    double s = TMath::Power(x[0] * 1E3, 2);

    const int nr = 4;
    double CR[nr] = {par[0], par[3], par[6], par[8]}; //1-par[0]-par[3]-par[6] };
    double CO[nr] = {par[1], par[4], par[7], par[8]}; //1-par[1]-par[4]-par[7] };
    double CP[nr] = {par[2], par[5], par[9]};         //1-par[2]-par[5], 0 };

    double KR[nr], KO[nr], KP[nr];
    for (int i = 0; i < nr; i++)
    {
        KR[i] = mode ? CR[i] / 2. : -CR[i] / 2.;
        KO[i] = CO[i] / 6.;
        KP[i] = CP[i] / 3.; //здесь уже нет дополнительного параметра n, как в F0 (в соотв. с моделью)
    }

    TComplex F1 = F0(x, par, mode);
    F1 += KR[1] * BW_Rho1(s) + KR[2] * BW_Rho2(s) + KR[3] * BW_Rho5(s);
    F1 += KO[1] * BW_Omg1(s) + KO[2] * BW_Omg2(s) + KO[3] * BW_Rho5(s);
    F1 += KP[1] * BW_Phi1(s) + KP[2] * BW_Phi2(s);

    return F1;
}
double MDVM::Cross_Section(double *x, double *par, bool mode)
{
    double s = TMath::Power(x[0] * 1E3, 2);
    if (x[0] < 0.4976 * 2)
        return 0;
    double fabs = TComplex::Abs(F1(x, par, mode));
    double M_K = mode ? m_kc : m_k0;
    double cs = (TMath::Pi() / 3.) * TMath::Power(ALPHA, 2) * C * TMath::Power(BETA(s, M_K), 3) * TMath::Power(fabs, 2) / s;
    return cs;
}

TF1 *MDVM::Cross_Section(bool mode)
{
    const int Npars = 11;
    TF1 *fcs_c = new TF1("Cross section", (mode ? Cross_Section_Charged : Cross_Section_Neutral), 0.98, 2.1, Npars);
    fcs_c->SetParNames("C_{#rho}", "C_{#omega}", "C_{#phi}", "C_{#rho(1450)}", "C_{#omega(1420)}", "C_{#phi(1680)}", "C_{#rho(1570)}", "C_{#omega(1650)}", "C_{#rho(2100)}", "C_{#phi(2180)}", "#eta");
    fcs_c->SetParameters(1.19, 1.42, 1, -0.092, -0.04, -0.114, -0.032, -0.105);
    return fcs_c;
}
