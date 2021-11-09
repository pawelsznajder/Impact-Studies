#include "HepMC3/GenEvent.h"

#include <TLorentzVector.h>

using namespace HepMC3;

Double_t scattered_Phi(TLorentzVector mup) {

	Double_t phi = mup.Phi();

	return phi;
}

Double_t scattered_Theta(TLorentzVector mup) {

        Double_t theta = mup.Theta();

        return theta;
}


double Phi_Phis(int mode, TLorentzVector q, TLorentzVector p, TLorentzVector mu, TLorentzVector mup, TLorentzVector v){

         //target spin
         TVector3 spin_tar;
         spin_tar.SetXYZ(0., 1., 0.);

         //variables
         TVector3 boost;
         TLorentzVector q_boosted, p_boosted, mu_boosted, mup_boosted, v_boosted, lvp1b2, lvp2b2;
         double sinb2, cosb2;
         double sign;

         //Phi
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

         //Phis
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

double Phi(TLorentzVector q, TLorentzVector p, TLorentzVector mu, TLorentzVector mup, TLorentzVector v){
        return Phi_Phis(0, q, p, mu, mup, v);
}

double Phis(TLorentzVector q, TLorentzVector p, TLorentzVector mu, TLorentzVector mup, TLorentzVector v){
        return Phi_Phis(1, q, p, mu, mup, v);
}

TLorentzVector getEIn(const GenEvent& evt){

        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){
                if((*it)->status() == 4 && (*it)->pid() == 11)
                        return TLorentzVector((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz(), (*it)->momentum().e());
        }

        std::cout << __func__ << ": error: no particle" << std::endl;
        exit(0);
}

TLorentzVector getEOut(const GenEvent& evt){

        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){
                if((*it)->status() == 1 && (*it)->pid() == 11 && (*it)->children().size() == 0)
                        return TLorentzVector((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz(), (*it)->momentum().e());
        }

        std::cout << __func__ << ": error: no particle" << std::endl;
        exit(0);
}

TLorentzVector getGammaStar(const GenEvent& evt){

        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){
                if((*it)->status() == 3 && (*it)->pid() == 22)
                        return TLorentzVector((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz(), (*it)->momentum().e());
        }

        std::cout << __func__ << ": error: no particle" << std::endl;
        exit(0);
}

TLorentzVector getPIn(const GenEvent& evt){

        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){
                if((*it)->status() == 4 && (*it)->pid() == 2212)
                        return TLorentzVector((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz(), (*it)->momentum().e());
        }

        std::cout << __func__ << ": error: no particle" << std::endl;
        exit(0);
}

TLorentzVector getPOut(const GenEvent& evt){

        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){
                if((*it)->status() == 1 && (*it)->pid() == 2212)
                        return TLorentzVector((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz(), (*it)->momentum().e());
        }

        std::cout << __func__ << ": error: no particle" << std::endl;
        exit(0);
}

TLorentzVector getGammaOut(const GenEvent& evt){

        const std::vector<ConstGenParticlePtr>& particles = evt.particles();

        for(std::vector<ConstGenParticlePtr>::const_iterator it = particles.begin(); it != particles.end(); it++){
                if((*it)->status() == 1 && (*it)->pid() == 22) {

                        for(std::vector<ConstGenParticlePtr>::const_iterator itV = (*it)->production_vertex()->particles_in().begin(); itV != (*it)->production_vertex()->particles_in().end(); itV++){
                                if((*itV)->status() == 4 && (*itV)->pid() == 2212)
                                        return TLorentzVector((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz(), (*it)->momentum().e());
                        }
                }
        }

        std::cout << __func__ << ": error: no particle" << std::endl;
        exit(0);
}

TLorentzVector getFourMomentum(std::shared_ptr<const HepMC3::GenParticle> p){
        return TLorentzVector(p->momentum().px(), p->momentum().py(), p->momentum().pz(), p->momentum().e());
}

class DVCSEvent{

        public:

        DVCSEvent(const GenEvent& evt, size_t mode){

                //particles
                TLorentzVector eIn;
                TLorentzVector eOut;
                TLorentzVector gammaStar;
                TLorentzVector pIn;
                TLorentzVector gammaOut;
                TLorentzVector pOut;

if(mode == 0){

                eIn = getFourMomentum(evt.particles().at(0));
                eOut = getFourMomentum(evt.particles().at(1));
                gammaStar = getFourMomentum(evt.particles().at(2));
                pIn = getFourMomentum(evt.particles().at(3));
                gammaOut = getFourMomentum(evt.particles().at(4));
                pOut = getFourMomentum(evt.particles().at(5));

}else{
                eIn = getEIn(evt);
                eOut = getEOut(evt);
                gammaStar = getGammaStar(evt);
		//gammaStar = getEIn(evt) - getEOut(evt);
                pIn = getPIn(evt);
                gammaOut = getGammaOut(evt);
                pOut = getPOut(evt);
}

                //variables

                m_Q2 = -1 * gammaStar.Mag2();

                m_xB = m_Q2 / (2 * pIn * gammaStar);

                m_y = (pIn * gammaStar) / (pIn * eIn);

                m_t = (pOut - pIn).Mag2();

                m_phi = Phi(gammaStar, pIn, eIn, eOut, gammaOut);

		m_eOut_Px = eOut.Px();

		m_eOut_Py = eOut.Py();

		m_eOut_Pz = eOut.Pz();

		m_eOut_E = eOut.E();

		m_gammaOut_Px = gammaOut.Px();

		m_gammaOut_Py = gammaOut.Py();

		m_gammaOut_Pz = gammaOut.Pz();

		m_gammaOut_E = gammaOut.E();

                scattered_phi = scattered_Phi(eOut);

		scattered_theta = scattered_Theta(eOut);
        }

        double getXb() const{
                return m_xB;
        }

        double getT() const{
                return m_t;
        }

        double getQ2() const{
                return m_Q2;
        }

        double getY() const{
                return m_y;
        }

        double getPhi() const{
                return m_phi;
        }

	double getEOut_Px() const{
		return m_eOut_Px;
	}

	double getEOut_Py() const{
                return m_eOut_Py;
        }

	double getEOut_Pz() const{
                return m_eOut_Pz;
        }

	double getEOut_E() const{
                return m_eOut_E;
        }

	double getGammaOut_Px() const{
                return m_gammaOut_Px;
        }

	double getGammaOut_Py() const{
                return m_gammaOut_Py;
        }

	double getGammaOut_Pz() const{
                return m_gammaOut_Pz;
        }

	double getGammaOut_E() const{
                return m_gammaOut_E;
        }

	Double_t getScatteredPhi() const{
		return scattered_phi;
	}

	Double_t getScatteredTheta() const{
		return scattered_theta;
	}

        private:

        double m_xB;
        double m_t;
        double m_Q2;
        double m_y;
        double m_phi;
	double m_eOut_Px;
	double m_eOut_Py;
	double m_eOut_Pz;
	double m_eOut_E;
	double m_gammaOut_Px;
	double m_gammaOut_Py;
	double m_gammaOut_Pz;
	double m_gammaOut_E;
	Double_t scattered_phi;
	Double_t scattered_theta;

};
