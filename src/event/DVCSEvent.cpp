#include "../../include/event/DVCSEvent.h"

#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <iostream>
#include <utility>

using namespace HepMC3;

DVCSEvent::DVCSEvent(const GenEvent& evt, int beamPolarisation, int beamCharge,
                const TVector3& targetPolarisation, bool isRCSample, int subProcessTypeMask) : BaseObject("DVCSEvent"){

        //four-momenta
        m_eIn = getEIn(evt);
        m_eOut = getEOut(evt);
        m_pIn = getPIn(evt);
        m_pOut = getPOut(evt);
        m_gammaOut = getGammaOut(evt);
        m_gammaStar = getGammaStar();

        m_gammaRC = getGammaRC(evt);

        //radiation
        m_rcTypeMask = 0;

        if(m_gammaRC.find(RCType::ISR) != m_gammaRC.end())
                m_rcTypeMask |= RCType::ISR;
        if(m_gammaRC.find(RCType::FSR) != m_gammaRC.end())
                m_rcTypeMask |= RCType::FSR;

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

void DVCSEvent::setFourMomentum(std::pair<TLorentzVector, TLorentzVector>& inputPair, KinematicsType::Type type, const TLorentzVector& mom) const{

        switch(type){

                case KinematicsType::Observed:{
                        inputPair.first = mom;
                        break;
                }

                case KinematicsType::True:{
                        inputPair.second = mom;
                        break;
                }

                default:{

                        std::cout << "error: " << __func__ << ": wrong kinematics type, " << 
                                KinematicsType(type).toString() << std::endl;
                        exit(0);
                }
        }
}

const TLorentzVector& DVCSEvent::getFourMomentum(const std::pair<TLorentzVector, TLorentzVector>& inputPair, KinematicsType::Type type) const{

        switch(type){

                case KinematicsType::Observed:{
                        return inputPair.first;
                        break;
                }

                case KinematicsType::True:{
                        return inputPair.second;
                        break;
                }

                default:{

                        std::cout << "error: " << __func__ << ": wrong kinematics type, " << 
                                KinematicsType(type).toString() << std::endl;
                        exit(0);
                }
        }
}

void DVCSEvent::printError(const std::string& functionName) const{

        std::cout << "error: " << functionName << ": inconsistent event" << std::endl;
        exit(0);
}

TLorentzVector DVCSEvent::makeTLorentzVector(const ConstGenParticlePtr& particle) const{
        return TLorentzVector(
                particle->momentum().px(), 
                particle->momentum().py(), 
                particle->momentum().pz(), 
                particle->momentum().e());
}

std::pair<TLorentzVector, TLorentzVector> DVCSEvent::getEIn(const GenEvent& evt) const{

        //to check if processes correctly
        bool isOK = false;

        //result
        std::pair<TLorentzVector, TLorentzVector> result;

        //check if have radiation
        bool hasRadiation = false;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only electron
                if((*it)->pid() != 11) continue;

                //ISR (is marked as beam and ends in e -> e + gamma vertex) 
                //TODO !!!
                if((*it)->status() == 4 && (*it)->children().size() == 1){

                        setFourMomentum(result, KinematicsType::Observed, makeTLorentzVector(*it));
                        hasRadiation = true;

                        // bool hasE = false;
                        // bool hasGamma = false;

                        // //reference and loop
                        // const std::vector<ConstGenParticlePtr>& particlesOut = (*it)->children();

                        // for(std::vector<ConstGenParticlePtr>::const_iterator itOut = particlesOut.begin(); itOut != particlesOut.end(); itOut++){

                        //         if((*itOut)->pid() == 11) hasE = true;
                        //         if((*itOut)->status() == 1 && (*itOut)->pid() == 22) hasGamma = true;
                        // } 

                        // //save
                        // if(hasE && hasGamma){
                        //         result.insert(std::make_pair(RCType::ISR, makeTLorentzVector(*it)));
                        //         hasRadiation = true;
                        // }
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
                        if(hasE && hasGammaStar){

                                setFourMomentum(result, KinematicsType::True, makeTLorentzVector(*it));
                                isOK = true;
                        }
                }
        }

        if(! hasRadiation){
                setFourMomentum(result, KinematicsType::Observed,  getFourMomentum(result, KinematicsType::True));   
        }

        if(! isOK) printError(__func__);
        return result;
}

std::pair<TLorentzVector, TLorentzVector> DVCSEvent::getEOut(const GenEvent& evt) const{

        //to check if processes correctly
        bool isOK = false;

        //result
        std::pair<TLorentzVector, TLorentzVector> result;

        //check if have radiation
        bool hasRadiation = false;

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
                                if(hasE && hasGamma){

                                        setFourMomentum(result, KinematicsType::Observed, makeTLorentzVector(*it));
                                        hasRadiation = true;
                                }
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
                                if(hasE && hasGammaStar){

                                        setFourMomentum(result, KinematicsType::True, makeTLorentzVector(*it));
                                        isOK = true;
                                }
                        }
                }
        }

        if(! hasRadiation){
                setFourMomentum(result, KinematicsType::Observed, getFourMomentum(result, KinematicsType::True));   
        }

        if(! isOK) printError(__func__);
        return result;
}

std::pair<TLorentzVector, TLorentzVector> DVCSEvent::getPIn(const GenEvent& evt) const{

        //to check if processes correctly
        bool isOK = false;

        //result
        std::pair<TLorentzVector, TLorentzVector> result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only proton
                if((*it)->pid() != 2212) continue;

                //Born (is beam) 
                if((*it)->status() == 4){

                        //save
                        setFourMomentum(result, KinematicsType::True, makeTLorentzVector(*it));
                        isOK = true;

                        setFourMomentum(result, KinematicsType::Observed, getFourMomentum(result, KinematicsType::True));
                }         
        }

        return result;
}

std::pair<TLorentzVector, TLorentzVector> DVCSEvent::getPOut(const GenEvent& evt) const{

        //to check if processes correctly
        bool isOK = false;

        //result
        std::pair<TLorentzVector, TLorentzVector> result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only proton
                if((*it)->pid() != 2212) continue;

                //Born (is not beam) 
                if((*it)->status() == 1){

                        //save
                        setFourMomentum(result, KinematicsType::True, makeTLorentzVector(*it));
                        isOK = true;

                        setFourMomentum(result, KinematicsType::Observed, getFourMomentum(result, KinematicsType::True));
                }         
        }

        if(! isOK) printError(__func__);
        return result;
}

std::pair<TLorentzVector, TLorentzVector> DVCSEvent::getGammaOut(const GenEvent& evt) const{

        //to check if processes correctly
        bool isOK = false;

        //result
        std::pair<TLorentzVector, TLorentzVector> result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only photon
                if((*it)->pid() != 22) continue;

                //Born (comes from gamma* + p -> gamma + p vertex) 
                if((*it)->parents().size() == 2){

                        bool hasP = false;
                        bool hasGammaStar = false;

                        //reference and loop
                        const std::vector<ConstGenParticlePtr>& particlesIn = (*it)->production_vertex()->particles_in();

                        for(std::vector<ConstGenParticlePtr>::const_iterator itIn = particlesIn.begin(); itIn != particlesIn.end(); itIn++){
                                
                                if((*itIn)->pid() == 2212) hasP = true;
                                if((*itIn)->status() == 13 && (*itIn)->pid() == 22) hasGammaStar = true;
                        }   

                        //save
                        if(hasP && hasGammaStar){

                                //save
                                setFourMomentum(result, KinematicsType::True, makeTLorentzVector(*it));
                                isOK = true;

                                setFourMomentum(result, KinematicsType::Observed, getFourMomentum(result, KinematicsType::True)); 
                        }   
                }        
        }

        if(! isOK) printError(__func__);
        return result;
}

std::pair<TLorentzVector, TLorentzVector> DVCSEvent::getGammaStar() const{

        //to check if processes correctly
        bool isOK = false;

        //result
        std::pair<TLorentzVector, TLorentzVector>  result;

        setFourMomentum(result, KinematicsType::True, getFourMomentum(m_eIn, KinematicsType::True) - getFourMomentum(m_eOut, KinematicsType::True));
        isOK = true;

        setFourMomentum(result, KinematicsType::Observed, getFourMomentum(m_eIn, KinematicsType::Observed) - getFourMomentum(m_eOut, KinematicsType::Observed));

        return result;
}

std::map<RCType::Type, TLorentzVector> DVCSEvent::getGammaRC(const GenEvent& evt) const{

         //to check if processes correctly
        bool isOK = true;

        //result
        std::map<RCType::Type, TLorentzVector> result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only photon
                if((*it)->pid() != 22) continue;

                //only real
                if((*it)->status() != 1) continue;

                //ISR
                if((*it)->parents().size() == 1){

                        if((*it)->parents().at(0)->pid() == 11 && (*it)->parents().at(0)->status() == 4){
                                 result.insert(std::make_pair(RCType::ISR, makeTLorentzVector(*it)));
                        }

                        //TODO
                        // if((*it)->parents().at(0)->pid() == 11 && (*it)->parents().at(0)->status() == 4 && (*it)->production_vertex()->particles_out().size() == 2){

                        //         bool hasE = false;
                        //         bool hasGamma = false;

                        //         //reference and loop
                        //         const std::vector<ConstGenParticlePtr>& particlesOut = (*it)->production_vertex()->particles_out();

                        //         for(std::vector<ConstGenParticlePtr>::const_iterator itOut = particlesOut.begin(); itOut != particlesOut.end(); itOut++){
                                        
                        //                 if((*itOut)->pid() == 11) hasE = true;
                        //                 if((*itOut)->pid() == 22) hasGamma = true;
                        //         }   

                        //         if(hasE && hasGamma){
                        //                 result.insert(std::make_pair(RCType::ISR, makeTLorentzVector(*it)));
                        //         }
                        // }
                }

                //FSR
                if((*it)->parents().size() == 1){

                    if((*it)->parents().at(0)->pid() == 11 && (*it)->parents().at(0)->status() == 1 && (*it)->production_vertex()->particles_out().size() == 2){
                           
                                bool hasE = false;
                                bool hasGamma = false;

                                //reference and loop
                                const std::vector<ConstGenParticlePtr>& particlesOut = (*it)->production_vertex()->particles_out();

                                for(std::vector<ConstGenParticlePtr>::const_iterator itOut = particlesOut.begin(); itOut != particlesOut.end(); itOut++){
                                        
                                        if((*itOut)->pid() == 11) hasE = true;
                                        if((*itOut)->pid() == 22) hasGamma = true;
                                }   


                                if(hasE && hasGamma){
                                        result.insert(std::make_pair(RCType::FSR, makeTLorentzVector(*it)));
                                }
                        }
                }
        }

        if(! isOK) printError(__func__);
        return result;
}

double DVCSEvent::getXB(KinematicsType::Type type) const{
        return -1 * getFourMomentum(m_gammaStar, type).Mag2() / (2 * getFourMomentum(m_pIn, type) * getFourMomentum(m_gammaStar, type));
}

double DVCSEvent::getT(KinematicsType::Type type) const{
        return (getFourMomentum(m_pOut, type) - getFourMomentum(m_pIn, type)).Mag2();
}

double DVCSEvent::getQ2(KinematicsType::Type type) const{
        return -1 * getFourMomentum(m_gammaStar, type).Mag2();
}

double DVCSEvent::getY(KinematicsType::Type type) const{
        return (getFourMomentum(m_pIn, type) * getFourMomentum(m_gammaStar, type)) / (getFourMomentum(m_pIn, type) * getFourMomentum(m_eIn, type));
}

double DVCSEvent::getPhi(KinematicsType::Type type) const{
        return getPhiPhiS(0, getFourMomentum(m_gammaStar, type), getFourMomentum(m_pIn, type), getFourMomentum(m_eIn, type), getFourMomentum(m_eOut, type), getFourMomentum(m_gammaOut, type));
}

double DVCSEvent::getPhiS(KinematicsType::Type type) const{
        return getPhiPhiS(1, getFourMomentum(m_gammaStar, type), getFourMomentum(m_pIn, type), getFourMomentum(m_eIn, type), getFourMomentum(m_eOut, type), getFourMomentum(m_gammaOut, type));
}

double DVCSEvent::getEtaEOut(KinematicsType::Type type) const{

        return getFourMomentum(m_eOut, type).Eta();
}

double DVCSEvent::getEtaPOut(KinematicsType::Type type) const{
        return getFourMomentum(m_pOut, type).Eta();
}

double DVCSEvent::getEtaGOut(KinematicsType::Type type) const{
        return getFourMomentum(m_gammaOut, type).Eta();
}

double DVCSEvent::getEGammaRC(RCType::Type type) const{

        std::map<RCType::Type, TLorentzVector>::const_iterator it = m_gammaRC.find(type);

        if(it == m_gammaRC.end()) return -1.E12;

        return it->second.E();
}

double DVCSEvent::getEtaGammaRC(RCType::Type type) const{

        std::map<RCType::Type, TLorentzVector>::const_iterator it = m_gammaRC.find(type);

        if(it == m_gammaRC.end()) return -1.E12;

        return it->second.Eta();
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

bool DVCSEvent::isRCSample() const{
        return m_isRCSample;
}

 bool DVCSEvent::checkRCType(int rcTypeMask) const{
        return (m_rcTypeMask & rcTypeMask);
 }

bool DVCSEvent::checkSubProcessType(int subProcessTypeMask) const{
        return (m_subProcessTypeMask & subProcessTypeMask);
}

        

  

