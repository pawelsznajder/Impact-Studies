#ifndef DVCSEVENT_H
#define DVCSEVENT_H

#include <HepMC3/GenEvent.h>
#include <TLorentzVector.h>

#include "../../include/other/BaseObject.h"

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
        double getXB() const;

        //get t
        double getT() const;

        //get Q2
        double getQ2() const;

        //get y
        double getY() const;

        //get phi
        double getPhi() const;

        //get phiS
        double getPhiS() const;
        
	//get pseudo-rapidity of outgoing electron
        double getEtaEOut() const;

	//get pseudo-rapidity of outgoing proton
        double getEtaPOut() const;

	//get pseudo-rapidity of outgoing photon
        double getEtaGOut() const;

        //get beam polarisation
        int getBeamPolarisation() const;

        //get beam charge
        int getBeamCharge() const;

        //get target polarisation
        TVector3 getTargetPolarisation() const;

        //get four-momentum of beam electron
        const TLorentzVector& getEIn() const;

        //get four-momentum of outgoing electron
        const TLorentzVector& getEOut() const;

        //get four-momentum of virtual photon
        const TLorentzVector& getGammaStar() const;

        //get four-momentum of beam proton
        const TLorentzVector& getPIn() const;

        //get four-momentum of outgoing proton
        const TLorentzVector& getPOut() const;

        //get four-momentum of real (DVCS) photon
        const TLorentzVector& getGammaOut() const;

private:

        //get phi (mode = 0) or phiS (mode = 1) angle
        double getPhiPhiS(int mode, const TLorentzVector& q, 
                const TLorentzVector& p, const TLorentzVector& mu, 
                const TLorentzVector& mup, const TLorentzVector& v) const;

        //find four-momentum of beam electron in HepMC3 event 
        TLorentzVector getEIn(const GenEvent& evt) const;

        //find four-momentum of outgoing electron in HepMC3 event 
        TLorentzVector getEOut(const GenEvent& evt) const;

        //find four-momentum of beam proton in HepMC3 event 
        TLorentzVector getPIn(const GenEvent& evt) const;

        //find four-momentum of outgoing proton in HepMC3 event 
        TLorentzVector getPOut(const GenEvent& evt) const;

        //find four-momentum of real (DVCS) photon in HepMC3 event 
        TLorentzVector getGammaOut(const GenEvent& evt) const;

        //four-momenta
        TLorentzVector m_eIn;
        TLorentzVector m_eOut;
        TLorentzVector m_gammaStar;
        TLorentzVector m_pIn;
        TLorentzVector m_pOut;
        TLorentzVector m_gammaOut;

        //variables
        double m_xB;
        double m_t;
        double m_Q2;
        double m_y;
        double m_phi;
        double m_phiS;
        double m_etaeOut;
        double m_etapOut;
        double m_etagOut;

        //beam polarisation
        int m_beamPolarisation;

        //beam charge
        int m_beamCharge;

        //target polarisation
        TVector3 m_targetPolarisation;
};

#endif
