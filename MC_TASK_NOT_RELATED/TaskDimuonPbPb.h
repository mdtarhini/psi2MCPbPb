/*
 *simple DIMuon Analysis task
 *
 */
#ifndef ALIANALYSISTASKJPSIPP_H
#define ALIANALYSISTASKJPSIPP_H

class TH1F;
class TH2F;
class TH3F;
class TNtuple;
class TTree;
class TList;
class TObjArray;
class TArrayF;
class TClonesArray;

class AliAnalysisManager;
class AliVEvent;
class AliMuonTrackCuts;
class AliAnalysisUtils;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class TaskDimuonPbPb : public AliAnalysisTaskSE {
public:
    TaskDimuonPbPb();
    TaskDimuonPbPb(const char *name);
    virtual ~TaskDimuonPbPb();
    virtual void     NotifyRun();
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
    Bool_t IsMuonFromJPsi(int muonIndex);

    private:

      AliMuonTrackCuts *fMuonTrackCuts;//!

      TClonesArray *fArrayMCParticles;//

      TTree *treeEvents; //!
      std::vector<double> vectorGenDimuon;
      std::vector<double> vectorRecDimuon;
      Int_t runNumber;
      Float_t centrality;

      TObjArray *fListEvent;

    ClassDef(TaskDimuonPbPb, 2);
};

#endif
