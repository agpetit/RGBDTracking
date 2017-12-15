#ifndef P_HELPER_INCLUDED
#define P_HELPER_INCLUDED

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h> 

#include <sys/time.h>

using namespace std;

namespace p_helper_socket{

	/*Called as:
		create_UDP_socket(sock, 6999, "192.168.1.2");

		Return Value:
			 0: No error
			-1: gethostbyname error (Unknown host)
			-2: connect function error

								     id socket, port,    ip of the server
	*/
	int socket_create_UDP_Sender(int &sock_id, int port, string ip_recipient){

		// Client Len
    	socklen_t clilen;
    	struct sockaddr_in server_str;

	    // Struct: Server, Client
	    struct sockaddr_in client_str; //server_str;

	    struct hostent *h;

    	server_str.sin_family=AF_INET;
    	server_str.sin_port=htons(port); // oppure atoi(port_unity_message.c_str()) se port fosse una stringa!
    	h=gethostbyname(ip_recipient.c_str());

	    if (h==0) {
	        return -1;
	    }

    	bcopy(h->h_addr,&server_str.sin_addr,h->h_length);
    	sock_id = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);

    	if(connect(sock_id, (struct sockaddr*) &server_str, sizeof(server_str))< 0){
    		return -2;
    	}

    	return 0;
	}



}



namespace p_helper_generic{

	// From String to Char
	inline char* string2char( string str) {	
		char *cstr = new char[str.length() + 1];
		strcpy(cstr, str.c_str());
		return cstr;
	}

	// converts a rate in Hz to an integer period in ms.
	inline uint16_t rateToPeriod(const double & rate) {
		if (rate > 0)
			return static_cast<uint16_t> (1000.0 / rate);
		else
			return 0;
	}
}



#endif
