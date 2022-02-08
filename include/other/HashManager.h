#ifndef HASHMANAGER_H
#define HASHMANAGER_H

#include <string>
#include <sstream>

/**
 * Singleton class that returnes unique strings of characters.
 */
class HashManager {

    public:

    	//get instance
        static HashManager* getInstance(){

            static HashManager instance;
            return &instance;
        }

        //get hash
        std::string getHash(){

			//get string of hash
			std::stringstream ss; 
			ss << (m_i++);

			//return
			return ss.str();
        }

    private:

    	//constructor
        HashManager(){
        	m_i = 0;
        }

        //iterator
        size_t m_i;

};

#endif