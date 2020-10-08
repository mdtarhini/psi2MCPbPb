/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

//include root libraries
#include <iostream>
#include <TCanvas.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TChain.h>
#include "TMath.h"
#include <TObjArray.h>
#include <TClonesArray.h>
//include aliroot libraries
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAODv0.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVZERO.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliAnalysisMuonUtility.h"
#include "AliAnalysisUtils.h"
#include "AliMuonTrackCuts.h"
#include "AliOADBMuonTrackCutsParam.h"
#include "AliCentrality.h"
#include "AliAODDimuon.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliAODHeader.h"
#include "AliESDMuonTrack.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
#include "AliAODDimuon.h"
#include "TaskDimuonPbPb.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"

ClassImp(TaskDimuonPbPb)

    //________________________________________________________________________
    TaskDimuonPbPb::TaskDimuonPbPb()
    : AliAnalysisTaskSE(),
      fArrayMCParticles(),
      fListEvent(0x0),
      fMuonTrackCuts(0),
      vectorGenDimuon(0),
      vectorRecDimuon(0),
      runNumber(0),
      centrality(300)
{
  // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
TaskDimuonPbPb::TaskDimuonPbPb(const char *name)
    : AliAnalysisTaskSE(name),
      fArrayMCParticles(),
      fListEvent(0x0),
      fMuonTrackCuts(0),
      vectorGenDimuon(0),
      vectorRecDimuon(0),
      runNumber(0),
      centrality(300)
{
  // Input slot #0 works with a TChain - it is connected to the default input container
  DefineInput(0, TChain::Class());
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TObjArray::Class()); // for output objarray
}

//________________________________________________________________________
TaskDimuonPbPb::~TaskDimuonPbPb()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fListEvent && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())
  {
    delete fListEvent;
  }
}

//________________________________________________________________________
void TaskDimuonPbPb::NotifyRun()
{
  /// Set run number for cuts
  if (fMuonTrackCuts)
    fMuonTrackCuts->SetRun(fInputHandler);
}

//________________________________________________________________________
void TaskDimuonPbPb::UserCreateOutputObjects()
{
  //Create the event, single muon, dimuon histgrams, to do that, it is only needed to add the name of the histogram to the corresponding list, it will be created and the bins will be set according to the name

  //Event histograms
  //Event histograms
  fListEvent = new TObjArray(2000);
  fListEvent->SetOwner(kTRUE);
  fListEvent->SetName("ListEvent");

  treeEvents = new TTree("eventsTree", "tree that contains information of the event");
  treeEvents->Branch("GenDimuon", &vectorGenDimuon);
  treeEvents->Branch("RecDimuon", &vectorRecDimuon);
  treeEvents->Branch("runNumber", &runNumber);
  treeEvents->Branch("centrality", &centrality);
  fListEvent->AddAtAndExpand(treeEvents, 0);

  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca | AliMuonTrackCuts::kMuMatchLpt);
  fMuonTrackCuts->SetIsMC();

  PostData(1, fListEvent);
}

//________________________________________________________________________
void TaskDimuonPbPb::UserExec(Option_t *)
{
  vectorGenDimuon.clear();
  vectorRecDimuon.clear();
  Float_t muonMass2 = AliAnalysisMuonUtility::MuonMass2();

  AliAODEvent *aod = 0;
  AliAODHeader *aodHeader = NULL;
  aod = static_cast<AliAODEvent *>(InputEvent());

  if (!aod)
  {
    AliError("ERROR: Could not retrieve ESD or AOD event !!");
    return;
  }

  TString strFiredTriggers; //MODIFICATIN FOR PbPb
  if (aod)
  {
    aodHeader = (AliAODHeader *)aod->GetHeader();
    strFiredTriggers = aod->GetFiredTriggerClasses();
  }
  else
  {
    AliError("ERROR: Could not retrieve ESD or AOD event !!");
    return;
  }

  UInt_t IsSelected =
      (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());

  if (!(strFiredTriggers.Contains("CINT7-B-NOPF-MUFAST") && (IsSelected & AliVEvent::kAny)))
  {
    return;
  }

  AliMultSelection *multSelection = (AliMultSelection *)aod->FindListObject("MultSelection");
  centrality = multSelection->GetMultiplicityPercentile("V0M", false);
  if (centrality > 90)
    return;

  runNumber = aod->GetRunNumber();

  fArrayMCParticles = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));

  for (int iParticle = 0; iParticle < fArrayMCParticles->GetEntries(); iParticle++)
  {

    AliAODMCParticle *particle = (AliAODMCParticle *)fArrayMCParticles->At(iParticle);

    if (particle->GetStatus() == 21)
      continue;

    if ((particle->GetPdgCode() == 100443))
    {
      vectorGenDimuon.push_back(particle->GetCalcMass());
      vectorGenDimuon.push_back(particle->Pt());
      vectorGenDimuon.push_back(particle->Y());
      break;
    }
  }

  //reconstructed
  TLorentzVector lvRecFirstMuon, lvRecSecondMuon, lvRecDimuon;
  int numberOfTracks = aod->GetNumberOfTracks();

  for (Int_t iFirstMuon = 0; iFirstMuon < numberOfTracks; iFirstMuon++)
  {

    AliVTrack *firstMuon = aod->GetTrack(iFirstMuon);

    if (!fMuonTrackCuts->IsSelected(firstMuon))
      continue;

    //test if the muon is coming from a decay of a Z-boson
    int firstMuonIndex = firstMuon->GetLabel();
    if (!IsMuonFromJPsi(firstMuonIndex))
      continue;
    for (Int_t iSecondMuon = iFirstMuon + 1; iSecondMuon < numberOfTracks; iSecondMuon++)
    {

      AliVTrack *secondMuon = aod->GetTrack(iSecondMuon);

      if (!fMuonTrackCuts->IsSelected(secondMuon))
        continue;

      int secondMuonIndex = secondMuon->GetLabel();
      if (!IsMuonFromJPsi(secondMuonIndex))
        continue;

      Float_t energy = muonMass2 + firstMuon->P() * firstMuon->P();
      energy = TMath::Sqrt(energy);
      lvRecFirstMuon.SetPxPyPzE(firstMuon->Px(), firstMuon->Py(), firstMuon->Pz(), energy);

      energy = muonMass2 + secondMuon->P() * secondMuon->P();
      energy = TMath::Sqrt(energy);
      lvRecSecondMuon.SetPxPyPzE(secondMuon->Px(), secondMuon->Py(), secondMuon->Pz(), energy);

      lvRecDimuon = lvRecFirstMuon + lvRecSecondMuon;

      if (firstMuon->Charge() == secondMuon->Charge())
        continue;

      vectorRecDimuon.push_back(lvRecDimuon.M());
      vectorRecDimuon.push_back(lvRecDimuon.Pt());
      vectorRecDimuon.push_back(lvRecDimuon.Rapidity());
    }
  }

  if (vectorGenDimuon.size() > 0)
    ((TTree *)fListEvent->UncheckedAt(0))->Fill();

  PostData(1, fListEvent);
}

//________________________________________________________________________
void TaskDimuonPbPb::Terminate(Option_t *)
{
  fListEvent = dynamic_cast<TObjArray *>(GetOutputData(1));
  if (!fListEvent)
  {
    AliError("Could not retrieve TObjArray* fListEvent");
    return;
  }
}

Bool_t TaskDimuonPbPb::IsMuonFromJPsi(int muonIndex)
{

  AliAODMCParticle *muonParticle = (AliAODMCParticle *)fArrayMCParticles->At(muonIndex);
  if (muonIndex < 0)
    return kFALSE;
  if ((muonParticle->GetPdgCode() == 13) || (muonParticle->GetPdgCode() == -13))
  {

    int iMother = muonParticle->GetMother();
    while (iMother >= 0)
    {
      AliAODMCParticle *firstMother = (AliAODMCParticle *)fArrayMCParticles->At(iMother);
      int pdgCodeOfFirstMother = firstMother->GetPdgCode();
      if (pdgCodeOfFirstMother == 100443)
      {
        return kTRUE;
      }
      iMother = firstMother->GetMother();
    }
  }
  return kFALSE;
}
