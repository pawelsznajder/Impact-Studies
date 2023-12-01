#include "../../include/event/DVCSEvent.h"

#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <iostream>
#include <utility>

using namespace HepMC3;

DVCSEvent::DVCSEvent(const GenEvent& evtRec, const GenEvent& evtGen, int beamPolarisation, int beamCharge,
                const TVector3& targetPolarisation, bool isRCSample, int subProcessTypeMask) : BaseObject("DVCSEvent"){

        //four-momenta
        for(size_t i = 0; i < 2; i++){

                KinematicsType::Type type;

                switch(i){
                        case 0: {type = KinematicsType::Observed; break;}
                        case 1: {type = KinematicsType::True; break;}
                }

                const GenEvent& evt = (i == 0)?(evtRec):(evtGen);

                setFourMomentum(m_eIn, type, getEIn(evt, type));
                setFourMomentum(m_eOut, type, getEOut(evt, type));
                setFourMomentum(m_pIn, type, getPIn(evt, type));
                setFourMomentum(m_pOut, type, getPOut(evt, type));
                setFourMomentum(m_gammaOut, type, getGammaOut(evt, type));

                setFourMomentum(m_gammaISR, type, getGammaISR(evt, type));
                setFourMomentum(m_gammaFSR, type, getGammaFSR(evt, type));
        }

        if(getFourMomentum(m_gammaISR, KinematicsType::True) != TLorentzVector()){
                setFourMomentum(m_eIn, KinematicsType::Born, getFourMomentum(m_eIn, KinematicsType::True) - getFourMomentum(m_gammaISR, KinematicsType::True));
        }else{
                setFourMomentum(m_eIn, KinematicsType::Born, getFourMomentum(m_eIn, KinematicsType::True));
        }

        if(getFourMomentum(m_gammaFSR, KinematicsType::True) != TLorentzVector()){
                setFourMomentum(m_eOut, KinematicsType::Born, getFourMomentum(m_eOut, KinematicsType::True) + getFourMomentum(m_gammaFSR, KinematicsType::True));
        }else{
                setFourMomentum(m_eOut, KinematicsType::Born, getFourMomentum(m_eOut, KinematicsType::True));
        }

        setFourMomentum(m_pIn, KinematicsType::Born, getFourMomentum(m_pIn, KinematicsType::True));
        setFourMomentum(m_pOut, KinematicsType::Born, getFourMomentum(m_pOut, KinematicsType::True));
        setFourMomentum(m_gammaOut, KinematicsType::Born, getFourMomentum(m_gammaOut, KinematicsType::True));

        setFourMomentum(m_gammaISR, KinematicsType::Born, TLorentzVector());
        setFourMomentum(m_gammaFSR, KinematicsType::Born, TLorentzVector());


        //in the case of reconstruction assign mass for charged particles and momentum for photons
        //check if reconstructed
        m_isReconstructed = true;

        const double c_electronMass = 0.510998910E-3;
        const double c_protonMass = 0.938272013;

        m_isReconstructed &= improveReconstruction(m_eIn, c_electronMass); 
        m_isReconstructed &= improveReconstruction(m_eOut, c_electronMass); 
        m_isReconstructed &= improveReconstruction(m_pIn, c_protonMass); 
        m_isReconstructed &= improveReconstruction(m_pOut, c_protonMass); 
        m_isReconstructed &= improveReconstruction(m_gammaOut, 0.); 

        improveReconstruction(m_gammaISR, 0.); 
        improveReconstruction(m_gammaFSR, 0.); 

        //four-momenta of virtual-photon
        for(size_t i = 0; i < 3; i++){

                KinematicsType::Type type;

                switch(i){
                        case 0: {type = KinematicsType::Observed; break;}
                        case 1: {type = KinematicsType::True; break;}
                        case 2: {type = KinematicsType::Born; break;}
                }

                setFourMomentum(m_gammaStar, type, getFourMomentum(m_eIn, type) - getFourMomentum(m_eOut, type));
        }

        //radiation
        m_rcTypeMask = 0;

        if(getFourMomentum(m_gammaISR, KinematicsType::True) != TLorentzVector())
                m_rcTypeMask |= RCType::ISR;
        if(getFourMomentum(m_gammaFSR, KinematicsType::True) != TLorentzVector())
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

void DVCSEvent::setFourMomentum(std::tuple<TLorentzVector, TLorentzVector, TLorentzVector>& input, KinematicsType::Type type, const TLorentzVector& mom) const{

        switch(type){

                case KinematicsType::Observed:{
                        std::get<0>(input) = mom;
                        break;
                }

                case KinematicsType::True:{
                        std::get<1>(input) = mom;
                        break;
                }

                case KinematicsType::Born:{
                        std::get<2>(input) = mom;
                        break;
                }

                default:{

                        std::cout << "error: " << __func__ << ": wrong kinematics type, " << 
                                KinematicsType(type).toString() << std::endl;
                        exit(0);
                }
        }
}

const TLorentzVector& DVCSEvent::getFourMomentum(const std::tuple<TLorentzVector, TLorentzVector, TLorentzVector>& input, KinematicsType::Type type) const{

        switch(type){

                case KinematicsType::Observed:{
                        return std::get<0>(input);
                        break;
                }

                case KinematicsType::True:{
                        return std::get<1>(input);
                        break;
                }

                case KinematicsType::Born:{
                        return std::get<2>(input);
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

TLorentzVector DVCSEvent::getEIn(const GenEvent& evt, KinematicsType::Type type) const{

        //to check if processes correctly
        bool isOK = false;

        //result
        TLorentzVector result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only electron
                if((*it)->pid() != 11) continue;

                //if observed just look for electron with status 4
                if((*it)->status() == 4){

                        result = makeTLorentzVector(*it);
                        isOK = true;
                        break;
                } 

                /*CHANGE
                //only electron
                if((*it)->pid() != 11) continue;

                //if observed just look for electron with status 4
                if(type == KinematicsType::Observed){
                        if((*it)->status() == 4){

                                result = makeTLorentzVector(*it);
                                isOK = true;
                                break;
                        } 
                }

                //if true look for electron with e -> e + gamma* vertex
                if(type == KinematicsType::True){

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

                                        result = makeTLorentzVector(*it);
                                        isOK = true;
                                        break;
                                }
                        }
                }
                */
        }

        //check if ok
        if(! isOK) printError(__func__);

        //return
        return result;
}

TLorentzVector DVCSEvent::getEOut(const GenEvent& evt, KinematicsType::Type type) const{

        //to check if processes correctly
        bool isOK = false;

        //result
        TLorentzVector result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only electron
                if((*it)->pid() != 11) continue;

                //if observed just look for electron with status 4
                if((*it)->status() == 1){

                        result = makeTLorentzVector(*it);
                        isOK = true;
                        break;
                } 

                /*CHANGE
                //only electron
                if((*it)->pid() != 11) continue;

                //if observed look for electron with no children
                if(type == KinematicsType::Observed){
                        if((*it)->children().size() == 0){

                                result = makeTLorentzVector(*it);
                                isOK = true;
                                break;
                        }
                }

                //if true look for electron that comes from e -> e + gamma* vertex 
                if(type == KinematicsType::True){
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

                                                result = makeTLorentzVector(*it);
                                                isOK = true;
                                                break;
                                        }
                                }
                        }
                }
                */
        }

        //check if ok
        if(! isOK) printError(__func__);

        //return
        return result;
}

TLorentzVector DVCSEvent::getPIn(const GenEvent& evt, KinematicsType::Type type) const{

        //to check if processes correctly
        bool isOK = false;

        //result
        TLorentzVector result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only proton
                if((*it)->pid() != 2212) continue;

                //is marked as beam 
                if((*it)->status() == 4){

                        result = makeTLorentzVector(*it);
                        isOK = true;
                        break;
                }  

                /*CHANGE
                //only proton
                if((*it)->pid() != 2212) continue;

                //is marked as beam 
                if((*it)->status() == 4){

                        result = makeTLorentzVector(*it);
                        isOK = true;
                        break;
                }
                */      
        }

        //check if ok
        if(! isOK) printError(__func__);

        //return
        return result;
}

TLorentzVector DVCSEvent::getPOut(const GenEvent& evt, KinematicsType::Type type) const{

        //to check if processes correctly
        bool isOK = false;

        //result
        TLorentzVector result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){


                //only proton
                if((*it)->pid() != 2212) continue;

                //is marked as beam 
                if((*it)->status() == 1){

                        result = makeTLorentzVector(*it);
                        isOK = true;
                        break;
                }  

                /*CHANGE
                //only proton
                if((*it)->pid() != 2212) continue;

                //is not marked as beam 
                if((*it)->status() == 1){

                        result = makeTLorentzVector(*it);
                        isOK = true;
                        break;
                }    
                */     
        }

        //check if ok
        if(! isOK) printError(__func__);

        //return
        return result;
}

TLorentzVector DVCSEvent::getGammaOut(const GenEvent& evt, KinematicsType::Type type) const{

        //to check if processes correctly
        bool isOK = false;

        //result
        TLorentzVector result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only photon
                if((*it)->pid() != 22) continue;

                //Comes from gamma* + p -> gamma + p vertex 
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

                                result = makeTLorentzVector(*it);
                                isOK = true;
                                break;
                        }   
                }        
        }

        //check if ok
        if(! isOK) printError(__func__);

        //return
        return result;
}

TLorentzVector DVCSEvent::getGammaISR(const GenEvent& evt, KinematicsType::Type type) const{

        //result
        TLorentzVector result;

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
                                result = makeTLorentzVector(*it);
                                break;
                        }
                }
        }

        //return
        return result;
}

TLorentzVector DVCSEvent::getGammaFSR(const GenEvent& evt, KinematicsType::Type type) const{

        //result
        TLorentzVector result;

        //reference and loop
        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){

                //only photon
                if((*it)->pid() != 22) continue;

                //only real
                if((*it)->status() != 1) continue;

                //FSR
                if((*it)->parents().size() == 1){

                           if((*it)->parents().at(0)->pid() == 11 && (*it)->parents().at(0)->status() == 3 && (*it)->production_vertex()->particles_out().size() == 2){
                           
                                bool hasE = false;
                                bool hasGamma = false;

                                //reference and loop
                                const std::vector<ConstGenParticlePtr>& particlesOut = (*it)->production_vertex()->particles_out();

                                for(std::vector<ConstGenParticlePtr>::const_iterator itOut = particlesOut.begin(); itOut != particlesOut.end(); itOut++){
                                        
                                        if((*itOut)->pid() == 11) hasE = true;
                                        if((*itOut)->pid() == 22) hasGamma = true;
                                }   

                                if(hasE && hasGamma){

                                        result = makeTLorentzVector(*it);
                                        break;
                                }
                        }
                }
        }

        //return
        return result;
}

const TLorentzVector& DVCSEvent::getEIn(KinematicsType::Type type) const{
        return getFourMomentum(m_eIn, type);
}

const TLorentzVector& DVCSEvent::getEOut(KinematicsType::Type type) const{
        return getFourMomentum(m_eOut, type);
}

const TLorentzVector& DVCSEvent::getPIn(KinematicsType::Type type) const{
        return getFourMomentum(m_pIn, type);
}

const TLorentzVector& DVCSEvent::getPOut(KinematicsType::Type type) const{
        return getFourMomentum(m_pOut, type);
}

const TLorentzVector& DVCSEvent::getGammaOut(KinematicsType::Type type) const{
        return getFourMomentum(m_gammaOut, type);
}

const TLorentzVector& DVCSEvent::getGammaISR(KinematicsType::Type type) const{

        if(type == KinematicsType::Born){

                std::cout << "error: " << __func__ << ": undefined for this type, " << KinematicsType(type).toString() << std::endl;
                exit(0);
        }

        return getFourMomentum(m_gammaISR, type);
}

const TLorentzVector& DVCSEvent::getGammaFSR(KinematicsType::Type type) const{

        if(type == KinematicsType::Born){
                
                std::cout << "error: " << __func__ << ": undefined for this type, " << KinematicsType(type).toString() << std::endl;
                exit(0);
        }

        return getFourMomentum(m_gammaFSR, type);
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

 bool DVCSEvent::checkRCType(RCType::Type rcType) const{
        return (m_rcTypeMask & rcType);
 }

bool DVCSEvent::checkSubProcessType(SubProcessType::Type subProcessType) const{
        return (m_subProcessTypeMask & subProcessType);
}

bool DVCSEvent::isReconstructed() const{
        return m_isReconstructed;
}    

bool DVCSEvent::improveReconstruction(std::tuple<TLorentzVector, TLorentzVector, TLorentzVector>& lvs, double mass) const{

        //check if reconstructed
        if(getFourMomentum(lvs, KinematicsType::Observed) == TLorentzVector()) return false;

        //get (TODO: copy made here)
        TLorentzVector lv = getFourMomentum(lvs, KinematicsType::Observed);

        //check if reconstructed
        if(lv.Px() == -1. && lv.Py() == -1. && lv.Pz() == -1. && lv.E() == -1.) return false;

        //assign energy
        lv.SetE(
                sqrt(
                        pow(mass, 2) + 
                        pow(lv.P(), 2)
                )
        );  
        

        //set
        setFourMomentum(lvs, KinematicsType::Observed, lv);

        //return
        return true;
}  
