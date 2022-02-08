#ifndef DVCSEVENT_H
#define DVCSEVENT_H

#include <HepMC3/GenEvent.h>
#include <TLorentzVector.h>

#include "../../include/other/BaseObject.h"
#include "../../include/other/RCType.h"

using namespace HepMC3;

/**
 * Reconstruction of DVCS event from HepMC3 event.
 */
class DVCSEvent : public BaseObject{

public:

        //constructor
        DVCSEvent(const GenEvent& evt, int beamPolarisation, int beamCharge,
                const TVector3& targetPolarisation);

        //destructor
        virtual ~DVCSEvent();

        //get xB
        double getXB(int rc = (RCType::ISR | RCType::FSR));

        //get t
        double getT(int rc = (RCType::ISR | RCType::FSR));

        //get Q2
        double getQ2(int rc = (RCType::ISR | RCType::FSR));

        //get y
        double getY(int rc = (RCType::ISR | RCType::FSR));

        //get phi
        double getPhi(int rc = (RCType::ISR | RCType::FSR));

        //get phiS
        double getPhiS(int rc = (RCType::ISR | RCType::FSR));
        
	//get pseudo-rapidity of outgoing electron
        double getEtaEOut(int rc = (RCType::ISR | RCType::FSR));

	//get pseudo-rapidity of outgoing proton
        double getEtaPOut(int rc = (RCType::ISR | RCType::FSR));

	//get pseudo-rapidity of outgoing photon
        double getEtaGOut(int rc = (RCType::ISR | RCType::FSR));

        //get beam polarisation
        int getBeamPolarisation() const;

        //get beam charge
        int getBeamCharge() const;

        //get target polarisation
        TVector3 getTargetPolarisation() const;

        //check if contains given type of radiation
        bool checkIfRC(RCType::Type rcType) const;

private:

        //get phi (mode = 0) or phiS (mode = 1) angle
        double getPhiPhiS(int mode, const TLorentzVector& q, 
                const TLorentzVector& p, const TLorentzVector& mu, 
                const TLorentzVector& mup, const TLorentzVector& v) const;

        //make TLorentzVector
        TLorentzVector makeTLorentzVector(const ConstGenParticlePtr& constGenParticlePtr) const;

        //find four-momentum of beam electron in HepMC3 event 
        std::map<RCType::Type, TLorentzVector> getEIn(const GenEvent& evt) const;

        //find four-momentum of outgoing electron in HepMC3 event 
        std::map<RCType::Type, TLorentzVector> getEOut(const GenEvent& evt) const;

        //find four-momentum of beam proton in HepMC3 event 
        std::map<RCType::Type, TLorentzVector> getPIn(const GenEvent& evt) const;

        //find four-momentum of outgoing proton in HepMC3 event 
        std::map<RCType::Type, TLorentzVector> getPOut(const GenEvent& evt) const;

        //find four-momentum of DVCS photon in HepMC3 event 
        std::map<RCType::Type, TLorentzVector> getGammaOut(const GenEvent& evt) const;

        //four-momenta
        std::map<RCType::Type, TLorentzVector> m_eIn;
        std::map<RCType::Type, TLorentzVector> m_eOut;
        std::map<RCType::Type, TLorentzVector> m_pIn;
        std::map<RCType::Type, TLorentzVector> m_pOut;
        std::map<RCType::Type, TLorentzVector> m_gammaOut;

        //this flag is used to avoid loading four-momenta for each evaluation of kinematic variables
        int m_this_rc;

        //load four-momenta for evaluation of kinematic variables
        void loadThisRC(int rc);

        //four-momenta for evaluation of kinematic variables
        TLorentzVector m_this_rc_eIn;
        TLorentzVector m_this_rc_eOut;
        TLorentzVector m_this_rc_pIn;
        TLorentzVector m_this_rc_pOut;
        TLorentzVector m_this_rc_gammaOut;

        TLorentzVector m_this_rc_gammaStar;

        //beam polarisation
        int m_beamPolarisation;

        //beam charge
        int m_beamCharge;

        //target polarisation
        TVector3 m_targetPolarisation;
};

#endif
