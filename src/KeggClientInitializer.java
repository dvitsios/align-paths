import java.io.IOException;
import keggapi.KEGGLocator;
import keggapi.KEGGPortType;


class KeggClientInitializer {

	public static KEGGLocator  locator;
	public static KEGGPortType serv;


	public KeggClientInitializer() throws IOException{
		try {
			locator = new KEGGLocator();
			serv = locator.getKEGGPort();
            
            
		}
		catch (javax.xml.rpc.ServiceException e) {
			System.err.println("Error: " + e.getMessage());
		}
	}
}

















