#include "../../include/event/DVCSEvent.h"

#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <iostream>

using namespace HepMC3;

DVCSEvent::DVCSEvent(const GenEvent& evt, int beamPolarisation, int beamCharge,
                const TVector3& targetPolarisation) : BaseObject("DVCSEvent"){

        //four-momenta
        m_eIn = getEIn(evt);
        m_eOut = getEOut(evt);
        m_gammaStar = getEIn(evt) - getEOut(evt);
        m_pIn = getPIn(evt);
        m_pOut = getPOut(evt);
        m_gammaOut = getGammaOut(evt);

        //variables
        m_Q2 = -1 * m_gammaStar.Mag2();
        m_xB = m_Q2 / (2 * m_pIn * m_gammaStar);
        m_y = (m_pIn * m_gammaStar) / (m_pIn * m_eIn);
        m_t = (m_pOut - m_pIn).Mag2();
        m_phi = getPhiPhiS(0, m_gammaStar, m_pIn, m_eIn, m_eOut, m_gammaOut);  
        m_phiS = getPhiPhiS(1, m_gammaStar, m_pIn, m_eIn, m_eOut, m_gammaOut);

        //beam polarisation
        m_beamPolarisation = beamPolarisation;

        //beam charge
        m_beamCharge = beamCharge;

        //target polarisation
        m_targetPolarisation = targetPolarisation;
}

DVCSEvent::~DVCSEvent(){
}

double DVCSEvent::getPhiPhiS(int mode, const TLorentzVector& q, 
                const TLorentzVector& p, const TLorentzVector& mu, 
                const TLorentzVector& mup, const TLorentzVector& v) const{

         //target spin
         TVector3 spin_tar;
         spin_tar.SetXYZ(0., 1., 0.);

         //variables
         TVector3 boost;
         TLorentzVector q_boosted, p_boosted, mu_boosted, mup_boosted, v_boosted, lvp1b2, lvp2b2;
         double sinb2, cosb2;
         double sign;

         //phi
         double Phi;
         boost = (q+p).BoostVector();
         q_boosted = q;                q_boosted.Boost(-boost);
         p_boosted = p;                p_boosted.Boost(-boost);
         mu_boosted = mu;              mu_boosted.Boost(-boost);
         mup_boosted = mup;            mup_boosted.Boost(-boost);
         v_boosted = v;                v_boosted.Boost(-boost);

         sign = ( ((mu_boosted.Vect()).Cross(mup_boosted.Vect())).Unit() ).Cross( ((q_boosted.Vect()).Cross(v_boosted.Vect())).Unit() ).Dot( (q_boosted.Vect()).Unit() );
         sign /= fabs(sign);

         sinb2 = ( ( ((mu_boosted.Vect()).Cross(mup_boosted.Vect())).Unit() ).Cross( ((q_boosted.Vect()).Cross(v_boosted.Vect())).Unit() ) ).Mag();
         cosb2 = ( ((mu_boosted.Vect()).Cross(mup_boosted.Vect())).Unit() ).Dot( ((q_boosted.Vect()).Cross(v_boosted.Vect())).Unit() );

         Phi = atan2(sign*sinb2, cosb2);
         if( Phi < 0. ) Phi = Phi + 2.*TMath::Pi();

         if( mode == 0 ) return Phi;

         //phiS
         double Phis;

         sign = ( ((mu.Vect()).Cross(mup.Vect())).Unit() ).Cross( (q.Vect()).Cross(spin_tar).Unit() ).Dot( (q.Vect()).Unit() );
         sign /= fabs(sign);

         sinb2 = ( ( ((mu.Vect()).Cross(mup.Vect())).Unit() ).Cross( (q.Vect()).Cross(spin_tar).Unit() ) ).Mag();
         cosb2 = ( ((mu.Vect()).Cross(mup.Vect())).Unit() ).Dot( (q.Vect()).Cross(spin_tar).Unit() );

         Phis = atan2(sign*sinb2, cosb2);
         if( Phis < 0. ) Phis = Phis + 2.*TMath::Pi();

         if( mode == 1 ) return Phis;

         return -1;
}

TLorentzVector DVCSEvent::getEIn(const GenEvent& evt) const{

        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){
                if((*it)->status() == 4 && (*it)->pid() == 11)
                        return TLorentzVector((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz(), (*it)->momentum().e());
        }

        std::cout << getClassName() << "::" << __func__ << 
                ": error: no particle" << std::endl;
        exit(0);
}

TLorentzVector DVCSEvent::getEOut(const GenEvent& evt) const{

        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){
                if((*it)->status() == 1 && (*it)->pid() == 11 && (*it)->children().size() == 0)
                        return TLorentzVector((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz(), (*it)->momentum().e());
        }

        std::cout << getClassName() << "::" << __func__ << 
                ": error: no particle" << std::endl;
        exit(0);
}

TLorentzVector DVCSEvent::getPIn(const GenEvent& evt) const{

        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){
                if((*it)->status() == 4 && (*it)->pid() == 2212)
                        return TLorentzVector((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz(), (*it)->momentum().e());
        }

        std::cout << getClassName() << "::" << __func__ << 
                ": error: no particle" << std::endl;
        exit(0);
}

TLorentzVector DVCSEvent::getPOut(const GenEvent& evt) const{

        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){
                if((*it)->status() == 1 && (*it)->pid() == 2212)
                        return TLorentzVector((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz(), (*it)->momentum().e());
        }

        std::cout << getClassName() << "::" << __func__ << 
                ": error: no particle" << std::endl;
        exit(0);
}

TLorentzVector DVCSEvent::getGammaOut(const GenEvent& evt) const{

        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){
                if((*it)->status() == 1 && (*it)->pid() == 22) {

                        for(std::vector<ConstGenParticlePtr>::const_iterator itV = (*it)->production_vertex()->particles_in().begin(); itV != (*it)->production_vertex()->particles_in().end(); itV++){
                                if((*itV)->status() == 4 && (*itV)->pid() == 2212)
                                        return TLorentzVector((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz(), (*it)->momentum().e());
                        }
                }
        }

        std::cout << getClassName() << "::" << __func__ << 
                ": error: no particle" << std::endl;
        exit(0);
}

double DVCSEvent::getXB() const{
        return m_xB;
}

double DVCSEvent::getT() const{
        return m_t;
}

double DVCSEvent::getQ2() const{
        return m_Q2;
}

double DVCSEvent::getY() const{
        return m_y;
}

double DVCSEvent::getPhi() const{
        return m_phi;
}

double DVCSEvent::getPhiS() const{
        return m_phiS;
}

const TLorentzVector& DVCSEvent::getEIn() const{
        return m_eIn;
}

const TLorentzVector& DVCSEvent::getEOut() const{
        return m_eOut;
}

const TLorentzVector& DVCSEvent::getGammaStar() const{
        return m_gammaStar;
}

const TLorentzVector& DVCSEvent::getPIn() const{
        return m_pIn;
}

const TLorentzVector& DVCSEvent::getPOut() const{
        return m_pOut;
}

const TLorentzVector& DVCSEvent::getGammaOut() const{
        return m_gammaOut;
}

int DVCSEvent::getBeamPolarisation() const{
        return m_beamPolarisation;
}

int DVCSEvent::getBeamCharge() const{
        return m_beamCharge;
}

TVector3 DVCSEvent::getTargetPolarisation() const{
        return m_targetPolarisation;
}