#ifndef DVCSEVENT_H
#define DVCSEVENT_H

#include <HepMC3/GenEvent.h>
#include <TLorentzVector.h>

#include "../../include/other/BaseObject.h"
#include "../../include/other/RCType.h"
#include "../../include/other/KinematicsType.h"


using namespace HepMC3;

/**
 * Reconstruction of DVCS event from HepMC3 event.
 */
class DVCSEvent : public BaseObject{

public:

        //constructor
        DVCSEvent(const GenEvent& evt, int beamPolarisation, int beamCharge,
                const TVector3& targetPolarisation, bool isRCSample, int subProcessTypeMask);

        //destructor
        virtual ~DVCSEvent();

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
        
	//get pseudo-rapidity of outgoing electron
        double getEtaEOut(KinematicsType::Type type = KinematicsType::Observed) const;

	//get pseudo-rapidity of outgoing proton
        double getEtaPOut(KinematicsType::Type type = KinematicsType::Observed) const;

	//get pseudo-rapidity of outgoing photon
        double getEtaGOut(KinematicsType::Type type = KinematicsType::Observed) const;

        //get energy of radiative photon of given type
        double getEGammaRC(RCType::Type type) const;

        //get pseudo-rapidity of radiative photon of given type
        double getEtaGammaRC(RCType::Type type) const;

        //get beam polarisation
        int getBeamPolarisation() const;

        //get beam charge
        int getBeamCharge() const;

        //get target polarisation
        TVector3 getTargetPolarisation() const;

        //check if from RC sample
        bool isRCSample() const;

        //check if contains given type of radiation
        bool checkRCType(int rcTypeMask) const;

        //check subprocess
        bool checkSubProcessType(int subProcessTypeMask) const;

private:

        //get phi (mode = 0) or phiS (mode = 1) angle
        double getPhiPhiS(int mode, const TLorentzVector& q, 
                const TLorentzVector& p, const TLorentzVector& mu, 
                const TLorentzVector& mup, const TLorentzVector& v) const;

        //make TLorentzVector
        TLorentzVector makeTLorentzVector(const ConstGenParticlePtr& constGenParticlePtr) const;

        //find four-momentum of beam electron in HepMC3 event 
        std::pair<TLorentzVector, TLorentzVector> getEIn(const GenEvent& evt) const;

        //find four-momentum of outgoing electron in HepMC3 event 
        std::pair<TLorentzVector, TLorentzVector> getEOut(const GenEvent& evt) const;

        //find four-momentum of beam proton in HepMC3 event 
        std::pair<TLorentzVector, TLorentzVector> getPIn(const GenEvent& evt) const;

        //find four-momentum of outgoing proton in HepMC3 event 
        std::pair<TLorentzVector, TLorentzVector> getPOut(const GenEvent& evt) const;

        //find four-momentum of DVCS photon in HepMC3 event 
        std::pair<TLorentzVector, TLorentzVector> getGammaOut(const GenEvent& evt) const;

        //evaluate four-momentum of virtual photon
        std::pair<TLorentzVector, TLorentzVector> getGammaStar() const;

        //find four-momenta of RC photons in HepMC3 event 
        std::map<RCType::Type, TLorentzVector> getGammaRC(const GenEvent& evt) const;

        //set four-momentum in a pair for specific kinematic type
        void setFourMomentum(std::pair<TLorentzVector, TLorentzVector>& inputPair, KinematicsType::Type type, const TLorentzVector& mom) const;

        //get reference of four-momentum corresponding to specific kinematic type
        const TLorentzVector& getFourMomentum(const std::pair<TLorentzVector, TLorentzVector>& inputPair, KinematicsType::Type type) const;

        //print error message and terminate program if processing of events is not successful 
        void printError(const std::string& functionName) const;
   
        //four-momenta
        std::pair<TLorentzVector, TLorentzVector> m_eIn;
        std::pair<TLorentzVector, TLorentzVector> m_eOut;
        std::pair<TLorentzVector, TLorentzVector> m_pIn;
        std::pair<TLorentzVector, TLorentzVector> m_pOut;
        std::pair<TLorentzVector, TLorentzVector> m_gammaOut;
        std::pair<TLorentzVector, TLorentzVector> m_gammaStar;

        std::map<RCType::Type, TLorentzVector> m_gammaRC;

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
};

#endif
