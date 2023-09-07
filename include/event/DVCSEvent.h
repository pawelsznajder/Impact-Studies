#ifndef DVCSEVENT_H
#define DVCSEVENT_H

#include <HepMC3/GenEvent.h>
#include <TLorentzVector.h>

#include "../../include/other/BaseObject.h"
#include "../../include/other/RCType.h"
#include "../../include/other/SubProcessType.h"
#include "../../include/other/KinematicsType.h"


using namespace HepMC3;

/**
 * Reconstruction of DVCS event from HepMC3 event.
 */
class DVCSEvent : public BaseObject{

public:

        //constructor
        DVCSEvent(const GenEvent& evtGen, const GenEvent& evtRec, int beamPolarisation, int beamCharge,
                const TVector3& targetPolarisation, bool isRCSample, int subProcessTypeMask);

        //destructor
        virtual ~DVCSEvent();

        //get four-momentum of beam electron  
        const TLorentzVector& getEIn(KinematicsType::Type type = KinematicsType::Observed) const;

        //get four-momentum of outgoing electron  
        const TLorentzVector& getEOut(KinematicsType::Type type = KinematicsType::Observed) const;

        //get four-momentum of beam proton 
        const TLorentzVector& getPIn(KinematicsType::Type type = KinematicsType::Observed) const;

        //get four-momentum of outgoing proton  
        const TLorentzVector& getPOut(KinematicsType::Type type = KinematicsType::Observed) const;

        //get four-momentum of DVCS photon 
        const TLorentzVector& getGammaOut(KinematicsType::Type type = KinematicsType::Observed) const;

        //get four-momenta of ISR photons  
        const TLorentzVector& getGammaISR(KinematicsType::Type type = KinematicsType::Observed) const;

        //get four-momenta of FSR photons  
        const TLorentzVector& getGammaFSR(KinematicsType::Type type = KinematicsType::Observed) const;

        //get xB
        double getXB(KinematicsType::Type type = KinematicsType::Observed) const;

        //get t
        double getT(KinematicsType::Type type = KinematicsType::Observed) const;

        //get Q2
        double getQ2(KinematicsType::Type type = KinematicsType::Observed) const;

        //get y
        double getY(KinematicsType::Type type = KinematicsType::Observed) const;

        //get phi
        double getPhi(KinematicsType::Type type = KinematicsType::Observed) const;

        //get phiS
        double getPhiS(KinematicsType::Type type = KinematicsType::Observed) const;

        //get beam polarisation
        int getBeamPolarisation() const;

        //get beam charge
        int getBeamCharge() const;

        //get target polarisation
        TVector3 getTargetPolarisation() const;

        //check if from RC sample
        bool isRCSample() const;

        //check if contains given type of radiation
        bool checkRCType(RCType::Type rcType) const;

        //check subprocess
        bool checkSubProcessType(SubProcessType::Type subProcessType) const;

        //check if reconstructed
        bool isReconstructed() const;

private:

        //get phi (mode = 0) or phiS (mode = 1) angle
        double getPhiPhiS(int mode, const TLorentzVector& q, 
                const TLorentzVector& p, const TLorentzVector& mu, 
                const TLorentzVector& mup, const TLorentzVector& v) const;

        //make TLorentzVector
        TLorentzVector makeTLorentzVector(const ConstGenParticlePtr& constGenParticlePtr) const;

        //find four-momentum of beam electron in HepMC3 event 
        TLorentzVector getEIn(const GenEvent& evt, KinematicsType::Type type) const;

        //find four-momentum of outgoing electron in HepMC3 event 
        TLorentzVector getEOut(const GenEvent& evt, KinematicsType::Type type) const;

        //find four-momentum of beam proton in HepMC3 event 
        TLorentzVector getPIn(const GenEvent& evt, KinematicsType::Type type) const;

        //find four-momentum of outgoing proton in HepMC3 event 
        TLorentzVector getPOut(const GenEvent& evt, KinematicsType::Type type) const;

        //find four-momentum of DVCS photon in HepMC3 event 
        TLorentzVector getGammaOut(const GenEvent& evt, KinematicsType::Type type) const;

        //find four-momenta of ISR photons in HepMC3 event 
        TLorentzVector getGammaISR(const GenEvent& evt, KinematicsType::Type type) const;

        //find four-momenta of FSR photons in HepMC3 event 
        TLorentzVector getGammaFSR(const GenEvent& evt, KinematicsType::Type type) const;

        //set four-momentum in a pair for specific kinematic type
        void setFourMomentum(std::pair<TLorentzVector, TLorentzVector>& inputPair, KinematicsType::Type type, const TLorentzVector& mom) const;

        //get reference of four-momentum corresponding to specific kinematic type
        const TLorentzVector& getFourMomentum(const std::pair<TLorentzVector, TLorentzVector>& inputPair, KinematicsType::Type type) const;

        //print error message and terminate program if processing of events is not successful 
        void printError(const std::string& functionName) const;

        //check if reconstructed and improve information
        bool improveReconstruction(std::pair<TLorentzVector, TLorentzVector>& lvs, double mass) const;
   
        //four-momenta
        std::pair<TLorentzVector, TLorentzVector> m_eIn;
        std::pair<TLorentzVector, TLorentzVector> m_eOut;
        std::pair<TLorentzVector, TLorentzVector> m_pIn;
        std::pair<TLorentzVector, TLorentzVector> m_pOut;
        std::pair<TLorentzVector, TLorentzVector> m_gammaOut;

        std::pair<TLorentzVector, TLorentzVector> m_gammaStar;

        std::pair<TLorentzVector, TLorentzVector> m_gammaISR;
        std::pair<TLorentzVector, TLorentzVector> m_gammaFSR;

        //beam polarisation
        int m_beamPolarisation;

        //beam charge
        int m_beamCharge;

        //target polarisation
        TVector3 m_targetPolarisation;

        //true if from sample including RCs
        bool m_isRCSample;

        //rc mask (indicated which types of RCs one observes in this event)
        int m_rcTypeMask;

        //subrpocess mask
        int m_subProcessTypeMask;

        //true is all basic particles are reconstructed (i.e. event reconstructed passes topology selections)
        bool m_isReconstructed;
};

#endif
