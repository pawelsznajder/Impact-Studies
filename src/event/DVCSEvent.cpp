#include "../../include/event/DVCSEvent.h"

#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <iostream>

using namespace HepMC3;

DVCSEvent::DVCSEvent(const GenEvent& evt, int beamPolarisation, int beamCharge,
                const TVector3& targetPolarisation, bool isRCSample, int subProcessTypeMask) : BaseObject("DVCSEvent"){

        //four-momenta
        m_eIn = getEIn(evt);
        m_eOut = getEOut(evt);
        m_pIn = getPIn(evt);
        m_pOut = getPOut(evt);
        m_gammaOut = getGammaOut(evt);

        //mask
        m_rcTypeMask = RCType::Born;

        if(m_eIn.find(RCType::ISR) != m_eIn.end()) m_rcTypeMask |= RCType::ISR;
        if(m_eOut.find(RCType::FSR) != m_eOut.end()) m_rcTypeMask |= RCType::FSR;

        //rc flag
        m_this_rc = -1;

        //beam polarisation
        m_beamPolarisation = beamPolarisation;

        //beam charge
        m_beamCharge = beamCharge;

        //target polarisation
        m_targetPolarisation = targetPolarisation;

        //from sample including RCs
        m_isRCSample = isRCSample;

        //sub process type
        m_subProcessTypeMask = subProcessTypeMask;
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


TLorentzVector DVCSEvent::makeTLorentzVector(const ConstGenParticlePtr& particle) const{
        return TLorentzVector(
                particle->momentum().px(), 
                particle->momentum().py(), 
                particle->momentum().pz(), 
                particle->momentum().e());
}

std::map<RCType::Type, TLorentzVector> DVCSEvent::getEIn(const GenEvent& evt) const{

        //result
        std::map<RCType::Type, TLorentzVector> result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only electron
                if((*it)->pid() != 11) continue;

                //ISR (is marked as beam and ends in e -> e + gamma vertex) 
                //TODO !!!
                if((*it)->status() == 4 && (*it)->children().size() == 1){

                        result.insert(std::make_pair(RCType::ISR, makeTLorentzVector(*it)));

                        // bool hasE = false;
                        // bool hasGamma = false;

                        // //reference and loop
                        // const std::vector<ConstGenParticlePtr>& particlesOut = (*it)->children();

                        // for(std::vector<ConstGenParticlePtr>::const_iterator itOut = particlesOut.begin(); itOut != particlesOut.end(); itOut++){

                        //         if((*itOut)->pid() == 11) hasE = true;
                        //         if((*itOut)->status() == 1 && (*itOut)->pid() == 22) hasGamma = true;
                        // } 

                        // //save
                        // if(hasE && hasGamma) result.insert(std::make_pair(RCType::ISR, makeTLorentzVector(*it)));
                }

                //Born (ends in e -> e + gamma* vertex) 
                if((*it)->children().size() == 2){

                        bool hasE = false;
                        bool hasGammaStar = false;

                        //reference and loop
                        const std::vector<ConstGenParticlePtr>& particlesOut = (*it)->children();

                        for(std::vector<ConstGenParticlePtr>::const_iterator itOut = particlesOut.begin(); itOut != particlesOut.end(); itOut++){

                                if((*itOut)->pid() == 11) hasE = true;
                                if((*itOut)->status() == 13 && (*itOut)->pid() == 22) hasGammaStar = true;
                        }   

                        //save
                        if(hasE && hasGammaStar) result.insert(std::make_pair(RCType::Born, makeTLorentzVector(*it)));
                }
        }

        return result;
}

std::map<RCType::Type, TLorentzVector> DVCSEvent::getEOut(const GenEvent& evt) const{

        //result
        std::map<RCType::Type, TLorentzVector> result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only electron
                if((*it)->pid() != 11) continue;

                //FSR (comes from e -> e + gamma vertex and does not end in any vertex) 
                if((*it)->parents().size() == 1 && (*it)->children().size() == 0){

                        if((*it)->production_vertex()->particles_out().size() == 2){

                                bool hasE = false;
                                bool hasGamma = false;

                                if((*it)->parents().at(0)->pid() == 11) hasE = true;

                                //reference and loop
                                const std::vector<ConstGenParticlePtr>& particlesOut = (*it)->production_vertex()->particles_out();

                                for(std::vector<ConstGenParticlePtr>::const_iterator itOut = particlesOut.begin(); itOut != particlesOut.end(); itOut++){
                                        if((*itOut)->status() == 1 && (*itOut)->pid() == 22) hasGamma = true;
                                }   

                                //save
                                if(hasE && hasGamma) result.insert(std::make_pair(RCType::FSR, makeTLorentzVector(*it)));
                        }
                }

                //Born (comes from e -> e + gamma* vertex) 
                if((*it)->parents().size() == 1){

                        if((*it)->production_vertex()->particles_out().size() == 2){

                                bool hasE = false;
                                bool hasGammaStar = false;

                                if((*it)->parents().at(0)->pid() == 11) hasE = true;

                                //reference and loop
                                const std::vector<ConstGenParticlePtr>& particlesOut = (*it)->production_vertex()->particles_out();

                                for(std::vector<ConstGenParticlePtr>::const_iterator itOut = particlesOut.begin(); itOut != particlesOut.end(); itOut++){
                                        if((*itOut)->status() == 13 && (*itOut)->pid() == 22) hasGammaStar = true;
                                }   

                                //save
                                if(hasE && hasGammaStar) result.insert(std::make_pair(RCType::Born, makeTLorentzVector(*it)));
                        }
                }
        }

        return result;
}

std::map<RCType::Type, TLorentzVector> DVCSEvent::getPIn(const GenEvent& evt) const{

        //result
        std::map<RCType::Type, TLorentzVector> result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only proton
                if((*it)->pid() != 2212) continue;

                //Born (is beam) 
                if((*it)->status() == 4){

                        //save
                        result.insert(std::make_pair(RCType::Born, makeTLorentzVector(*it)));
                }         
        }

        return result;
}

std::map<RCType::Type, TLorentzVector> DVCSEvent::getPOut(const GenEvent& evt) const{

        //result
        std::map<RCType::Type, TLorentzVector> result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only proton
                if((*it)->pid() != 2212) continue;

                //Born (is not beam) 
                if((*it)->status() == 1){

                        //save
                        result.insert(std::make_pair(RCType::Born, makeTLorentzVector(*it)));
                }         
        }

        return result;
}

std::map<RCType::Type, TLorentzVector> DVCSEvent::getGammaOut(const GenEvent& evt) const{

        //result
        std::map<RCType::Type, TLorentzVector> result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only proton
                if((*it)->pid() != 22) continue;

                //Born (comes from gamma* + p -> gamma + p vertex) 
                if((*it)->parents().size() == 2){

                        bool hasP = false;
                        bool hasGammaStar = false;

                        //reference and loop
                        const std::vector<ConstGenParticlePtr>& particlesOut = (*it)->production_vertex()->particles_in();

                        for(std::vector<ConstGenParticlePtr>::const_iterator itOut = particlesOut.begin(); itOut != particlesOut.end(); itOut++){
                                
                                if((*itOut)->pid() == 2212) hasP = true;
                                if((*itOut)->status() == 13 && (*itOut)->pid() == 22) hasGammaStar = true;
                        }   

                        //save
                        if(hasP && hasGammaStar) result.insert(std::make_pair(RCType::Born, makeTLorentzVector(*it)));     
                }        
        }

        return result;
}

double DVCSEvent::getXB(int rc) {

        loadThisRC(rc);

        return -1 * m_this_rc_gammaStar.Mag2() / (2 * m_this_rc_pIn * m_this_rc_gammaStar);
}

double DVCSEvent::getT(int rc) {

        loadThisRC(rc);

        return (m_this_rc_pOut - m_this_rc_pIn).Mag2();
}

double DVCSEvent::getQ2(int rc) {

        loadThisRC(rc);

        return -1 * m_this_rc_gammaStar.Mag2();
}

double DVCSEvent::getY(int rc) {

        loadThisRC(rc);

        return (m_this_rc_pIn * m_this_rc_gammaStar) / (m_this_rc_pIn * m_this_rc_eIn);
}

double DVCSEvent::getPhi(int rc) {

        loadThisRC(rc);

        return getPhiPhiS(0, m_this_rc_gammaStar, m_this_rc_pIn, m_this_rc_eIn, m_this_rc_eOut, m_this_rc_gammaOut);
}

double DVCSEvent::getPhiS(int rc) {

        loadThisRC(rc);

        return getPhiPhiS(1, m_this_rc_gammaStar, m_this_rc_pIn, m_this_rc_eIn, m_this_rc_eOut, m_this_rc_gammaOut);
}

double DVCSEvent::getEtaEOut(int rc) {

        loadThisRC(rc);

        return m_this_rc_eOut.Eta();
}

double DVCSEvent::getEtaPOut(int rc) {

        loadThisRC(rc);

        return m_this_rc_pOut.Eta();
}

double DVCSEvent::getEtaGOut(int rc) {

        loadThisRC(rc);

        return m_this_rc_gammaOut.Eta();
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

void DVCSEvent::loadThisRC(int rc){

        //check if the same
        if(rc == m_this_rc) return;

        //iterator
        std::map<RCType::Type, TLorentzVector>::const_iterator it;

        //eIn
        it = m_eIn.end();

        if(rc & RCType::ISR){
                it = m_eIn.find(RCType::ISR);  
        }

        if(it == m_eIn.end()){
                it = m_eIn.find(RCType::Born);
        }

        m_this_rc_eIn = it->second;

        //eOut
        it = m_eOut.end();

        if(rc & RCType::FSR){
                it = m_eOut.find(RCType::FSR);  
        }

        if(it == m_eOut.end()){
                it = m_eOut.find(RCType::Born);
        }

        m_this_rc_eOut = it->second;

        //pIn
        it = m_pIn.find(RCType::Born);
        
        m_this_rc_pIn = it->second;

        //pOut
        it = m_pOut.find(RCType::Born);
        
        m_this_rc_pOut = it->second;

        //gamma DVCS
        it = m_gammaOut.find(RCType::Born);
        
        m_this_rc_gammaOut = it->second;

        //gamma *
        m_this_rc_gammaStar = m_this_rc_eIn - m_this_rc_eOut;

        //save
        m_this_rc = rc;
}

bool DVCSEvent::isRCSample() const{
        return m_isRCSample;
}

 bool DVCSEvent::checkRCType(int rcTypeMask) const{
        return (m_rcTypeMask & rcTypeMask);
 }

bool DVCSEvent::checkSubProcessType(int subProcessTypeMask) const{
        return (m_subProcessTypeMask & subProcessTypeMask);
}

        

  

