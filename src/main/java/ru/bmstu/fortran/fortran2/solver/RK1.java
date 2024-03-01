package ru.bmstu.fortran.fortran2.solver;

import static ru.bmstu.fortran.fortran2.solver.common.CGED.*;
import static ru.bmstu.fortran.fortran2.solver.common.CGEN.*;
import static ru.bmstu.fortran.fortran2.solver.common.CNUMBR.*;
import static ru.bmstu.fortran.fortran2.solver.common.CVLE.*;
import static ru.bmstu.fortran.fortran2.solver.common.CVLT.*;
import static ru.bmstu.fortran.fortran2.solver.common.CVNT.*;
import static ru.bmstu.fortran.fortran2.solver.common.CW.*;
import static ru.bmstu.fortran.fortran2.solver.common.XK.*;
import static ru.bmstu.fortran.fortran2.solver.common.XL.*;


public class RK1 {
    public static void main(String[] args) {
        trueInit();

    }

    private static void rk2() {
        double[][][] A = new double[46][22][16];
        int N1 = 46;
        int N2 = 22;
        int N3 = 16;
        INM = IN - 1;
        JNM = JN - 1;
        NITER = 0;
        NVAR = 0;
        do {
            NVAR++;
            ALFA += 0.2;
            double DELD = 0.5;
            if (NVAR != 1)
                DELD *= 1.E-3;
            grid(N1, N2, N3, A);
            C = 3.14 * DOOK * DELOK;
            double S1 = 0.785 * Math.pow(DCMB, 2);
            double G1 = 3.14 * DOGOR * DELGOR;
            for (int I = 1; I <= IN; I++) {
                PKS[I - 1] = P;
                if (I < IC1) {
                    continue;
                }
                int JL = JC + IC - I;
                double DSS = X2[JL - 1] * 2;
                double FSS = 3.14 / 4 * Math.pow(DSS, 2);
                double FSKR = 3.14 * Math.pow(X2[JC - 1], 2);
                double FOS = FSS / FSKR;
                if (FOS <= 1)
                    FOS = 1;
                double FA2 = 1 / FOS;
                double APR = Math.pow(2.0 / (AKT + 1.0), AKT / (AKT - 1.0)) * P;
                double AP2 = APR / P;
                double AP1 = 1;
                int NIT = 0;
                double Z1;
                double Z2;
                double APQ;
                double ZAP1;
                double ZAP2;
                double ZAP3;
                double ZAP4;
                double FA1;
                double FZ;
                // 82
                while (true) {
                    NIT++;
                    Z1 = AP1;
                    Z2 = AP2;
                    APQ = (Z1 + Z2) / 2;
                    ZAP1 = Math.pow(2.0 / (AKT + 1), (AKT + 1.0) / (2.0 * (AKT - 1.0)));
                    ZAP2 = Math.pow((AKT - 1) / 2, 0.5);
                    ZAP3 = Math.pow(1 / APQ, 1 / AKT);
                    ZAP4 = Math.pow(1.0 - Math.pow(APQ, (AKT - 1.0) / AKT), 0.5);
                    FA1 = 1 / (ZAP1 * ZAP2 * ZAP3 / ZAP4);
                    FZ = FA2 - FA1;
                    if (Math.abs(FZ) <= 1.E-6 || NIT > 1.E+4) {
                        break;
                    }
                    if (FA2 <= FA1)
                        AP2 = APQ;
                    if (FA2 >= FA1)
                        AP1 = APQ;
                }
                PKS[I - 1] = APQ * P;
                if (I > IC)
                    PKS[I - 1] = PKS[IC];
            }
            double APZ = P - PKS[IC1 - 1];
            for (int I = 2; I <= IC1; I++) {
                PKS[I - 1] = P - APZ * X1[I - 1] / X1[IC1 - 1];
                double APO = PKS[I - 1] / P;
            }
            double BETA = 0;
            if (ALFA <= 1.0) {
                BETA = 1794.0 + 491.252 * ALFA - 475.9317 * Math.pow(ALFA, 2);
            }

            if (ALFA > 1.0) {
                BETA = 2154.733 - 402.252 * ALFA + 54.1698 * Math.pow(ALFA, 2);
            }

            if (ALFA < 0.6) {
                BETA = 1917.0 - 2500.0 * (0.6 - ALFA);
            }

            if (BETA < 600.0) {
                BETA = 600.0;
            }
            double AKM1 = 70.0;
            DKR = 2.0 * X2[JC];
            double AKM = 4.0 * ALFA;
            double FKR = Math.pow(DKR, 2) * 3.14 / 4.0;
            AMS = P * FKR / BETA;
            AMG = AMS / (1 + AKM);
            AMO = AMG * AKM;
            double AMO2 = AMO;
            ROS = P / 77427.4;
            ROP = P / 154855.7;
            double AMS1 = D;
            double AKW = S1 / (C + S1);
            AMS1 = AKW * AMO;
            VINS1 = AMS1 / (ROS * 3.14 * Math.pow(R[JZ], 2));
            double AMO1 = AMS1;
            double AMG1 = 0.0;
            VINS = (AMO - AMO1) / ROS / (3.14 * (Math.pow(R[JA2], 2) - Math.pow(R[JA1], 2)));
            VINP = (AMG - AMG1) / ROP / (3.14 * (Math.pow(R[JA], 2) - Math.pow(R[JZ1], 2)));
            DKS = X2[JB1] * 2.0;
            ALKS = X1[IC1] - X1[IC2];
            ALKD = X1[IC];
            double ALPRIW = (Math.pow(R[JA3], 2) * X1[IC2] + (Math.pow(R[JB1], 2) - Math.pow(R[JC], 2)) * (X1[IC1] - X1[IC2]) + (Math.pow(R[JB1], 2) - Math.pow(R[JC], 2)) * (X1[IC] - X1[IC1]) / 3.0 + Math.pow(R[JC], 2) * (X1[IC] - X1[IC2])) / Math.pow(R[JC], 2);
            TGY = (X2[JN] - X2[JC]) / (X1[IC] - X1[IC1]);
            DKAR = X2[JA3] * 2.0;
            DELKAR = X1[IC2];
            DOOK = X2[JA1] * 2.0;
            DELOK = X2[JA2] - X2[JA1];
            DOGOR = X2[JZ1] * 2.0;
            DELGOR = X2[JA] - X2[JZ1];
            DCMB = X2[JZ] * 2.0;
            init(N1, N2, N3, A);
            double VR = A[1 - 1][JZ1 + 1 - 1][NMU2] / R[JZ1 + 1 - 1];

            String header = "7 " + String.format("%25s", "RASHET N" + NVAR);
            String line1 = String.format(" P=%11.4f ALFA=%6.2f KM=%6.2f BETA=%6.1f Ms=%11.4f Mok=%11.4f Mgor=%11.4f Vgor=%6.1f Vok=%6.1f ROg=%6.1f ROok=%6.1f Vokr=%6.1f Vcmw=%6.1f Mcmw=%11.4f",
                    P, ALFA, AKM, BETA, AMS, AMO, AMG, VINP, VINS, ROP, ROS, VR, VINS1, AMS1);
            String line2 = String.format(" JZ=%3d JZ1=%3d JA=%3d JA1=%3d JA2=%3d JA3=%3d JC=%3d JN=%3d IC=%3d IC1=%3d IC2=%3d IN=%3d",
                    JZ, JZ1, JA, JA1, JA2, JA3, JC, JN, IC, IC1, IC2, IN);
            String line3 = String.format(" DKS=%11.4f ALKS=%11.4f ALKD=%11.4f DKR=%11.4f tgY=%6.1f Dprk=%11.4f Lprk=%11.4f Dook=%11.4f DELOK=%11.4f DOGOR=%11.4f DELGOR=%11.4f Dcmw=%11.4f Lpriw=%11.4f Fod=%11.4f Fcmb=%11.4f Fgor=%11.4f Kwosp=%11.4f",
                    DKS, ALKS, ALKD, DKR, TGY, DKAR, DELKAR, DOOK, DELOK, DOGOR, DELGOR, DCMB, ALPRIW, C, S1, G1, AKW);
            System.out.println(header);
            System.out.println(line1);
            System.out.println(line2);
            System.out.println(line3);

            NITER = 0;

            do {
                NITER++;
                viscos(N1, N2, N3, A);
                veldis(N1, N2, N3, A);
                convec(N1, N2, N3, A);
                eqn(N1, N2, N3, A);
            } while (NITER < NMAX);

            System.out.println(String.format("%25sNITER=%4d", "", NITER));

            int JC17 = JC - 1;

            double ZGCP = 0.0;
            double ZHP = 0.0;
            if (ALFA > 1.0) {
                for (int J = 2; J <= JC17; J++) {
                    ZHP += A[IN - 1][J - 1][NMU1 - 1] * (A[IN - 1][J + 1 - 1][NF - 1] - A[IN - 1][J - 1 - 1][NF - 1]) / 2.0;
                    ZGCP += A[IN - 1][J - 1][NMFU1 - 1] * (A[IN - 1][J + 1 - 1][NF - 1] - A[IN - 1][J - 1 - 1][NF - 1]) / 2.0;
                }

                ZGCP += A[IN - 1][1 - 1][NMFU1 - 1] * (A[IN - 1][2 - 1][NF - 1] - A[IN - 1][1 - 1][NF - 1]) / 2.0
                        + A[IN - 1][JC - 1][NMFU1 - 1] * (A[IN - 1][JC - 1][NF - 1] - A[IN - 1][JC17 - 1][NF - 1]) / 2.0;

                ZGCP /= (A[1 - 1][JA - 1][NF - 1] - A[1 - 1][JZ1 - 1][NF - 1]) + 0.013833 * (A[1 - 1][JZ - 1][NF - 1] - A[1 - 1][1 - 1][NF - 1]);
                ZHP /= (A[1 - 1][JC - 1][NF - 1] - A[1 - 1][1 - 1][NF - 1]);

                if (Math.abs(ZGCP) > 1.0) {
                    ZGCP = 1.0;
                }

                double ZGCPE = Math.sqrt(1.0 - Math.abs(ZGCP));

                System.out.printf(" ************* KOEF Fkomp (GOR) = %.5f **********%n", ZGCPE);
            } else {
                ZGCP = 0;
                for (int J = 2; J <= JC17; J++) {
                    ZGCP = ZGCP + A[IN - 1][J - 1][NMOX1 - 1] * (A[IN - 1][J + 1 - 1][NF - 1] - A[IN - 1][J - 1 - 1][NF - 1]) / 2;
                }

                ZGCP += A[IN - 1][1 - 1][NMOX1 - 1] * (A[IN - 1][2 - 1][NF - 1] - A[IN - 1][1 - 1][NF - 1]) / 2.0
                        + A[IN - 1][JC - 1][NMOX1 - 1] * (A[IN - 1][JC - 1][NF - 1] - A[IN - 1][JC17 - 1][NF - 1]) / 2.0;

                ZGCP *= 6.28 / AKM;
                ZGCP /= AMG;

                if (Math.abs(ZGCP) > 1)
                    ZGCP = 1;
                double ZGCPT = Math.sqrt(1 - Math.abs(ZGCP));
                System.out.printf(" ************* KOEF Fkomp (OK) = %.5f **********%n", ZGCPT);
            }
            int IL = INM;

            double Z = 0.0;

            for (int J = 2; J <= JC; J++) {
                IL = INM;
                double ZO = A[IL - 1][J - 1 - 1][NMOX1 - 1] + 4.27 * A[IL - 1][J - 1 - 1][NMPR1 - 1] + 3.35 * A[IL - 1][J - 1 - 1][NMOX2 - 1] +
                        376.0 * A[IL - 1][J - 1 - 1][NMDF - 1] + 2.785 * A[IL - 1][J - 1 - 1][NMPR2 - 1] + 2.765 * A[IL - 1][J - 1 - 1][NMFU2 - 1];
                double ZG = A[IL - 1][J - 1 - 1][NMFU1 - 1] + 1.09 * A[IL - 1][J - 1 - 1][NMPR1 - 1] + 0.857 * A[IL - 1][J - 1 - 1][NMOX2 - 1] +
                        96.0 * A[IL - 1][J - 1 - 1][NMDF - 1] + 0.711 * A[IL - 1][J - 1 - 1][NMPR2 - 1] + 0.706 * A[IL - 1][J - 1 - 1][NMFU2 - 1];
                double ALFAJ = ZO / ZG / 3.9166;

                double BETAJ = 0;
                if (ALFAJ <= 1.0) {
                    BETAJ = 1794. + 491.252 * ALFAJ - 475.9317 * Math.pow(ALFAJ, 2);
                }
                if (ALFAJ > 1) {
                    BETAJ = 2154. - 402.252 * ALFAJ + 54.1698 * Math.pow(ALFAJ, 2);
                }

                if (ALFAJ < 0.6) {
                    BETAJ = 1917. - 2500. * (0.6 - ALFAJ);
                }

                if (BETAJ < 600.0) {
                    BETAJ = 600.0;
                }

                double ZF = (A[IL - 1][J - 1][NF - 1] - A[IL - 1][J - 1 - 1][NF - 1]) / (A[IL - 1][JC - 1][NF - 1] - A[IL - 1][1 - 1][NF - 1]);
                Z += BETAJ * ZF;

            } //RK 246
            double FICM = Z / BETA;
            print1(N1, N2, N3, A, IV);

            double FIOP = A[1 - 1][JA - 1][NF - 1] - A[1 - 1][JZ1 - 1][NF - 1] + 0.013833 * (A[1 - 1][JZ - 1][NF - 1] - A[1 - 1][1 - 1][NF - 1]);
            double FIOS = A[IC2 - 1][JA2 - 1][NF - 1] - A[IC2 - 1][JA1 - 1][NF - 1] + 0.986167 * (A[1 - 1][JZ - 1][NF - 1] - A[1 - 1][1 - 1][NF - 1]);

            FIOP = FIOP * 3.1416 * 2.0;
            FIOS = FIOS * 3.1416 * 2.0;
            double DOP = FIOS / FIOP;
        } while (NVAR != NVARM);
    }

    private static void print1(int N1, int N2, int N3, double[][][] A, int IV) {
    }

    private static void eqn(int N1, int N2, int N3, double[][][] A) {
        double SOURCE = 0;
        if (!(NITER <= NPRINT && NVAR == 1)) {
            RL = 0.1;
            for (int J = 2; J <= JNM; J++) {
                int IL = IMIN[J - 1];
                int IM = IMAX[J - 1];
                if (J == JA3)
                    IL = IC2 + 1;
                if (J == JC)
                    IM = IC - 1;
                for (int I = IL; I <= IM; I++) {
                    if ((J == JA && I == 2) ||
                            (J == JA1 && I == (IC2 + 1)) ||
                            (J == JA2 && I == (IC2 + 1)) ||
                            (J == JA3 && I == (IC2 + 1)) ||
                            (J == JZ && I == 2) ||
                            (J == JZ1 && I == 2)) {
                        //goto 111
                    }
                    if (JB1 != JN && J == (JB1 - 1) && I == IC2) {
                        //goto 113
                    }
                    double Z = A[I - 1][J - 1][NW - 1];
                    SOURCE = sorce(N1, N2, N3, A, SOURCE, I, J, NW);
                    double RSO = R[J - 1] * R[J - 1];
                    double BBE = 2. * RSO * BE[I - 1];
                    double BBW = 2. * RSO * BW[I - 1];
                    double BBS = (R[J - 1 - 1] * R[J - 1 - 1] + RSO) * BS[J - 1];
                    double BBN = (R[J + 1 - 1] * R[J + 1 - 1] + RSO) * BN[J - 1];
                    double AAE = RSO * AE[I - 1][J - 1];
                    double AAW = RSO * AW[I - 1][J - 1];
                    double AAN = RSO * AN[I - 1][J - 1];
                    double AAS = RSO * AS[I - 1][J - 1]; //RK3 44
                    double ANUM = (AAE + A[I + 1 - 1][J - 1][NMU - 1] * BBE) * A[I + 1 - 1][J - 1][NW - 1] +
                            (AAW + A[I - 1 - 1][J - 1][NMU - 1] * BBW) * A[I - 1 - 1][J - 1][NW - 1] +
                            (AAN + A[I - 1][J + 1 - 1][NMU - 1] * BBN) * A[I - 1][J + 1 - 1][NW - 1] +
                            (AAS + A[I - 1][J - 1 - 1][NMU - 1] * BBS) * A[I - 1][J - 1 - 1][NW - 1] +
                            SOURCE;
                    double ADNM = AAE + AAW + AAN + AAS + A[I - 1][J - 1][NMU - 1] * (BBE + BBW + BBN + BBS);
                    if (ADNM == 0) {
                        //goto 112
                    }
                    A[I - 1][J - 1][NW - 1] = ANUM / ADNM;
                    A[I - 1][J - 1][NW - 1] = Z + (A[I - 1][J - 1][NW - 1] - Z) * RL;
                    //112
                    //goto 114
                    //111
                    Z = A[I - 1][J - 1][NW - 1];
                    double RISO = 1 / R[J - 1] / R[J - 1];
                    double ROP3 = A[I - 1][J - 1][NRO3 - 1];
                    BBE = 4 / (A[I + 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BE[I - 1];
                    BBW = 4 / (A[I - 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BW[I - 1];
                    BBN = 16 / (A[I - 1][J + 1 - 1][NRO3 - 1] + ROP3) / Math.pow(R[J + 1 - 1] + R[J - 1], 2) * BN[J - 1];
                    BBS = 16 / (A[I - 1][J - 1 - 1][NRO3 - 1] + ROP3) / Math.pow(R[J - 1 - 1] + R[J - 1], 2) * BS[J - 1];
                    ANUM = BBE * A[I + 1 - 1][J - 1][NF - 1] + BBW * A[I - 1 - 1][J - 1][NF - 1] + BBN * A[I - 1][J + 1 - 1][NF - 1] + BBS * A[I - 1][J - 1 - 1][NF - 1];
                    ADNM = BBE + BBW + BBN + BBS;
                    A[I - 1][J - 1][NW - 1] = A[I - 1][J - 1][NF - 1] * ADNM - ANUM;
                    A[I - 1][J - 1][NW - 1] = Z + (A[I - 1][J - 1][NW - 1] - Z) * RL;
                    //goto 114
                    //113
                    Z = A[I - 1][J - 1][NW - 1];
                    RISO = 1 / R[J - 1] / R[J - 1];
                    ROP3 = A[I - 1][J - 1][NRO3 - 1];
                    BBE = 4 / (A[I + 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BE[I - 1];
                    BBW = 4 / (A[I - 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BW[I - 1];
                    BBN = 16 / (A[I - 1][J + 1 - 1][NRO3 - 1] + ROP3) / Math.pow(R[J + 1 - 1] + R[J - 1], 2) * BN[J - 1];
                    BBS = 16 / (A[I - 1][J - 1 - 1][NRO3 - 1] + ROP3) / Math.pow(R[J - 1 - 1] + R[J - 1], 2) * BS[J - 1];
                    ANUM = BBE * A[I + 1 - 1][J - 1][NF - 1] + BBW * A[I - 1 - 1][J - 1][NF - 1] + BBN * A[I - 1][J + 1 - 1][NF - 1] + BBS * A[I - 1][J - 1 - 1][NF - 1];
                    ADNM = BBE + BBW + BBN + BBS;
                    A[I - 1][J - 1][NW - 1] = A[I - 1][J - 1][NF - 1] * ADNM - ANUM;
                    A[I - 1][J - 1][NW - 1] = Z + (A[I - 1][J - 1][NW - 1] - Z) * RL;
                    //114
                    if (Math.abs(A[I - 1][J - 1][NW - 1]) >= 1.E-20) {
                        //goto 15
                    }
                    double RS = 0;
                    if (Math.abs(Z) >= 1.E-20) {
                        RS = 1 - A[I - 1][J - 1][NW - 1] / Z;
                    }
                    //goto 16
                    RS = 1 - Z / A[I - 1][J - 1][NW - 1]; // 15
                    //16
                    if (Math.abs(A[I - 1][J - 1][NW - 1]) < 1.E-20 && Math.abs(Z) < 1.E-20) {
                        RS = 0;
                    }
                    if (Math.abs(RS) > Math.abs(RSDU[NW - 1])) {
                        RSDU[NW - 1] = RS;
                    }
                    if (NITER >= 20) {
                        A[I - 1][J - 1][NW - 1] = Z + (A[I - 1][J - 1][NW - 1] - Z);
                    }
                    if (A[I - 1][J - 1][NW - 1] >= 1.E+14) {
                        A[I - 1][J - 1][NW - 1] = 1.E+14;
                    }
                    if (A[I - 1][J - 1][NW - 1] <= -1.E+14) {
                        A[I - 1][J - 1][NW - 1] = -1.E+14;
                    }
                }
            } //11
        } // 177 RK3 93
    }

    private static double sorce(int N1, int N2, int N3, double[][][] A, double SOURCE, int I, int J, int Z) {
        return 0;
    }

    private static void convec(int N1, int N2, int N3, double[][][] A) {

    }

    private static void veldis(int N1, int N2, int N3, double[][][] A) {

    }

    private static void viscos(int N1, int N2, int N3, double[][][] A) {

    }

    private static void grid(int N1, int N2, int N3, double[][][] A) {
    }

    private static void init(int N1, int N2, int N3, double[][][] A) {
        for (int J = 1; J <= JA3; J++) {
            IMIN[J - 1] = 2;
        }
        int JAA = JA3 + 1;
        for (int J = JAA; J <= JN; J++) {
            IMIN[J - 1] = IC2 + 1;
        }
        for (int J = 1; J <= JC; J++) {
            IMAX[J - 1] = IN - 1;
        }
        for (int J = JC1; J <= JN; J++) {
            IMAX[J - 1] = IC - 1 + JC - J;
        }
        for (int J = JA1; J <= JA2; J++) {
            A[IC2 - 1][J - 1][NV1 - 1] = VINS;
            A[IC2 - 1][J - 1][NMFU1 - 1] = 0;
            A[IC2 - 1][J - 1][NRO3 - 1] = ROS;
            A[IC2 - 1][J - 1][NMOX1 - 1] = ZMO2N;
            A[IC2 - 1][J - 1][NT - 1] = 298;
            A[IC2 - 1][J - 1][NMU1 - 1] = 0;
            if (NCORD == 2)
                A[IC2 - 1][J - 1][NF - 1] = ROS * VINS * (R[J - 1] - R[JA2 - 1]) * 0.5 * (R[J - 1] + R[JA2 - 1]);
            if (NCORD == 1)
                A[IC2 - 1][J - 1][NF - 1] = A[IC2 - 1][JA2 - 1][NF - 1] - ROS * VINS * (X2[JA2 - 1] - X2[J - 1]);
        }
        A[IC2 - 1][JA2 - 1][NV1 - 1] = 0;
        for (int J = JA3; J <= JA1; J++) {
            A[IC2 - 1][J - 1][NF - 1] = A[IC2 - 1][JA1 - 1][NF - 1];
        }
        for (int I = 1; I <= IC2; I++) {
            A[I - 1][JA3 - 1][NF - 1] = A[IC2 - 1][JA3 - 1][NF - 1];
        }
        A[IC2 - 1][JA1 - 1][NV1 - 1] = 0;
        for (int J = JA; J <= JA3; J++) {
            A[1 - 1][J - 1][NF - 1] = A[1 - 1][JA3 - 1][NF - 1];
        }
        for (int J = JZ1; J <= JA; J++) {
            A[1 - 1][J - 1][NV1 - 1] = VINP;
            A[1 - 1][J - 1][NMU2 - 1] = VINP * Math.pow(R[J - 1], 2) / R[JA - 1] * 1.5;
            A[1 - 1][J - 1][NMOX1 - 1] = 0;
            A[1 - 1][J - 1][NMFU1 - 1] = ZMF2N;
            A[1 - 1][J - 1][NRO3 - 1] = ROP;
            A[1 - 1][J - 1][NT - 1] = 298;
            A[1 - 1][J - 1][NMU1 - 1] = H;
            if (NCORD == 2)
                A[1 - 1][J - 1][NF - 1] = A[1 - 1][JA - 1][NF - 1] + ROP * VINP * (R[J - 1] - R[JA - 1]) * 0.5 * (R[J - 1] + R[JA - 1]);
            if (NCORD == 1)
                A[1 - 1][J - 1][NF - 1] = A[1 - 1][JA - 1][NF - 1] + ROP * VINP * (X2[J - 1] - X2[JA - 1]);
        }
        A[1 - 1][JA - 1][NV1 - 1] = 0;
        A[1 - 1][JZ1 - 1][NV1 - 1] = 0;
        for (int J = JZ; J <= JZ1; J++) {
            A[1 - 1][J - 1][NF - 1] = A[1 - 1][JZ1 - 1][NF - 1];
        }
        for (int J = 1; J <= JZ; J++) {
            A[1 - 1][J - 1][NT - 1] = 298;
            A[1 - 1][J - 1][NV1 - 1] = VINS1;
            A[1 - 1][J - 1][NMOX1 - 1] = 1;
            A[1 - 1][J - 1][NRO3 - 1] = 1. / (A[1 - 1][J - 1][NMFU1 - 1] / 16. + A[1 - 1][J - 1][NMOX1 - 1] / 32. + A[1 - 1][J - 1][NMPR1 - 1] / 44.
                    + A[1 - 1][J - 1][NMOX2 - 1] / 28. + A[1 - 1][J - 1][NMPR2 - 1] / 18. + A[1 - 1][J - 1][NMFU2 - 1] / 17.
                    + A[1 - 1][J - 1][NMDF - 1] / 2.);
            A[1 - 1][J - 1][NRO3 - 1] = A[1 - 1][J - 1][NRO3 - 1] * P / 8314.4 / A[1 - 1][J - 1][NT - 1];
            if (NCORD == 2)
                A[1 - 1][J - 1][NF - 1] = A[1 - 1][JZ - 1][NF - 1] + VINS1 * A[1 - 1][J - 1][NRO3 - 1] *
                        (R[J - 1] - R[JZ - 1]) * 0.5 * (R[J - 1] + R[JZ - 1]);
        }
        for (int I = 1; I <= IN; I++) {
            A[I - 1][1 - 1][NF - 1] = A[1 - 1][1 - 1][NF - 1];
        }
        A[1 - 1][JZ - 1][NV1 - 1] = 0;
        for (int J = JZ1; J <= JA; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1] + 1;
            if (J <= JC)
                IM = IL;
            if (J == JC1 || J == JB)
                IM = IMAX[J - 1] + 1;
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NMFU1 - 1] = 1;
                A[I - 1][J - 1][NMOX1 - 1] = 0;
                A[I - 1][J - 1][NRO3 - 1] = ROP;
                A[I - 1][J - 1][NT - 1] = 298;
                A[I - 1][J - 1][NMU1 - 1] = H;
            }
        }
        int JA17 = JA + 1;
        int JA18 = JA1 - 1;
        for (int J = JA17; J <= JA18; J++) {
            int IL = IMIN[J - 1] - 1;
            int IM = IMAX[J - 1] + 1;
            if (J <= JC)
                IM = IN;
            if (J == JC1 || J == JB)
                IM = IMAX[J - 1] + 1;
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NMOX1 - 1] = 1;
                A[I - 1][J - 1][NMFU2 - 1] = 0;
                A[I - 1][J - 1][NRO3 - 1] = ROS;
                A[I - 1][J - 1][NT - 1] = 298;
                A[I - 1][J - 1][NMFU1 - 1] = 0;
                A[I - 1][J - 1][NMU1 - 1] = H / 100;
            }
        }
        JA17 = JA1 + 1;
        for (int J = JA1; J <= JN; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1] + 1;
            if (J <= JC)
                IM = IN;
            if (J == JC1 || J == JB)
                IM = IMAX[J - 1] + 1;
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NMFU1 - 1] = 0;
                A[I - 1][J - 1][NT - 1] = 298;
                A[I - 1][J - 1][NRO3 - 1] = ROS;
                A[I - 1][J - 1][NMFU2 - 1] = 0;
                A[I - 1][J - 1][NMOX1 - 1] = 1;
                A[I - 1][J - 1][NMU1 - 1] = H / 100;
            }
        }
        for (int J = 1; J <= JZ; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1] + 1;
            if (J <= JC)
                IM = IN;
            if (J == JB || J == JC1)
                IM = IMAX[J - 1] + 1;
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NMOX1 - 1] = A[1 - 1][J - 1][NMOX1 - 1];
                A[I - 1][J - 1][NRO3 - 1] = A[1 - 1][J - 1][NRO3 - 1];
                A[I - 1][J - 1][NMU1 - 1] = A[1 - 1][J - 1][NMU1 - 1];
                A[I - 1][J - 1][NMPR1 - 1] = A[1 - 1][J - 1][NMPR1 - 1];
                A[I - 1][J - 1][NMPR2 - 1] = A[1 - 1][J - 1][NMPR2 - 1];
                A[I - 1][J - 1][NMFU1 - 1] = A[1 - 1][J - 1][NMFU1 - 1];
                A[I - 1][J - 1][NT - 1] = ROM1 + 10;
            }
        }
        int JZ17 = JZ1 - 1;
        int JZ18 = JZ + 1;
        for (int J = JZ18; J <= JZ17; J++) {
            int IL = IMIN[J - 1] - 1;
            int IM = IMAX[J - 1] + 1;
            if (J <= JC)
                IM = IN;
            if (J == JB || J == JC1)
                IM = IMAX[J - 1] + 1;
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NMOX1 - 1] = 0;
                A[I - 1][J - 1][NMFU1 - 1] = 1;
                A[I - 1][J - 1][NRO3 - 1] = ROP;
                A[I - 1][J - 1][NMU1 - 1] = H;
                A[I - 1][J - 1][NT - 1] = 298;
            }
        }
        int JAB = JA3 - 1;
        for (int J = 2; J <= JAB; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1];
            if (J == JC)
                IM = IC - 1;
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NF - 1] = A[1 - 1][J - 1][NF - 1];
            }
        } //RK1 215
        for (int J = JA3; J <= JNM; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1];
            if (J == JC)
                IM = IC - 1;
            if (J == JA3)
                IL = IC2 + 1;
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NF - 1] = A[IC2 - 1][J - 1][NF - 1];
            }
        }
    }

    private static void trueInit() {
        DKR = 10.0;
        DCMB = 2.0;
        DOGOR = 8.0;
        DELGOR = 1.0;
        DKAR = 12.0;
        DOOK = 30.0;
        DELOK = 1.5;
        DKS = 24.0;
        DELKAR = 10.0;
        ALKS = 36.0;
        ALKD = 42.239;
        TGY = 0.573;
        ALFA = 0.1;
        NW = 1;
        NF = 2;
        NMFU1 = 3;
        NMOX1 = 4;
        NMPR1 = 5;
        NMPR2 = 6;
        NMDF = 7;
        NMOX2 = 8;
        NMFU2 = 9;
        NMU1 = 10;
        NMU2 = 16;
        NRO3 = 12;
        NT = 13;
        NV1 = 14;
        NV2 = 15;
        NMU = 11;
        NK = 17;
        IE = 10;
        IV = 16;
        NMAX = 1000;
        NPRINT = 200;
        IP = 1;
        CC = 0.0035;
        NVARM = 1;
        VINP = 9.329155;
        VINS = 26.03549;
        ZMF2N = 1.00;
        ZMO2N = 1.00;
        VINT = 0.1;
        VINS1 = 2.E-3;
        ROP = 6.458;
        ROS = 12.916;
        RO1 = 6.458;
        ROM1 = 815.0;
        ROT = 0.1;
        NITER = 0;
        LN = 20;
        AFI = 0.34;
        AKT = 1.2;
        D = 5.E-3;
        P = 10.E+5;
        H = -4681.06;
        NCORD = 2;
        JA = 9;
        JA1 = 15;
        JC1 = 13;
        JB = 21;
        IC1 = 31;
        JA2 = 18;
        JB1 = 22;
        IC = 45;
        IN = 46;
        JN = 22;
        IC2 = 10;
        JZ = 4;
        JZ1 = 6;
        JA3 = 12;
        JC = 6;
        ROS1 = 6.458;
        ROP1 = 12.916;
        VINP1 = 0.1;
        VINT1 = 0.1;
        X1 = new double[]{0.0, 0.6, 1.2, 1.8, 2.4, 3.0, 4.4, 5.8, 7.2, 8.6, 10.0,
                11.4, 12.8, 14.2, 15.6, 17.0, 18.4, 19.8, 21.2, 22.6, 24.0,
                25.4, 26.8, 28.2, 29.6, 31.0, 32.4, 33.8, 35.2, 36.6, 38.0,
                39.4, 40.8, 42.2, 43.6, 45.0, 46.4, 47.8, 49.2, 50.6, 52.0,
                53.0, 54.0, 55.0, 56.0, 57.0};
        X2 = new double[]{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.6, 7.1, 8.0, 8.75, 9.5,
                10.25, 11.0, 12.0, 13.0, 14.0, 15.0, 15.5, 16.0, 16.5, 17.0};

        RSDU = new double[12];
        LMIN = new int[25];
    }
}
