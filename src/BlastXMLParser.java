/*
 *  Copyright (C) 2012 by Vitsios Dimitrios
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in
 *  all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */


import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.StringTokenizer;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;

class BlastXMLParser {

	//max expectation value accepted for succesful alignment filtering.
	public static Double expValueThreshold;
	
	public BlastXMLParser(){

	}

	public void parseXMLFile(String xmlFile, String baseOrg, String queryOrg, int baseIdx, int queryIdx){
		
		int tmpCnt = 1;
		int hitsCounter = 0;
		int baseOffset = 0;
		int queryOffset = 0;
		for (int i=0; i<queryIdx; i++) {
			queryOffset += MainClass.pr[i].length;
		}
		for (int i=0; i<baseIdx; i++) {
			baseOffset += MainClass.pr[i].length;
		}
		
		File fXmlFile = new File(xmlFile);
		DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
		try {
			DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
			Document doc = dBuilder.parse(fXmlFile);

			doc.getDocumentElement().normalize();
			
			//System.out.println("\n---------------------------------------------------");
			//System.out.println("Output from "+xmlFile+" parsing:\n");
							   
			NodeList iterationList = doc.getElementsByTagName("Iteration");
			//"s" is the queryGeneIdx.
			for (int s = 0; s < iterationList.getLength(); s++) {
				Node iterNode = iterationList.item(s);
				if (iterNode.getNodeType() == Node.ELEMENT_NODE) {
					Element iterElement = (Element) iterNode;
					//System.out.println("Query gene definition: "  + getTagValue("Iteration_query-def",eElement));
					NodeList iterHitsList = iterElement.getElementsByTagName("Iteration_hits");
					for (int t = 0; t < iterHitsList.getLength(); t++) {
						Node iterHitsNode = iterHitsList.item(t);
						if (iterHitsNode.getNodeType() == Node.ELEMENT_NODE) {
							Element iterHitsElement = (Element) iterHitsNode;
							NodeList hitList = iterHitsElement.getElementsByTagName("Hit");
							double minExpVal = 1000;
							for (int u = 0; u < hitList.getLength(); u++) {
								int baseGeneIdx = -1;
								Node hitNode = hitList.item(u);
								if (hitNode.getNodeType() == Node.ELEMENT_NODE) {
									Element hitElement = (Element) hitNode;
									//get hit_gene string
									String hitGeneStr = getTagValue("Hit_def",hitElement);
									
									StringTokenizer stok = new StringTokenizer(hitGeneStr);
									String hitGeneName = stok.nextToken();
									for (int fg=0; fg < MainClass.pr[baseIdx].length; fg++) {
										if(MainClass.pr[baseIdx][fg].equals(hitGeneName)){
											baseGeneIdx = fg;
											break;
										}
									}				
									//System.out.println("HitGeneString: "+hitGeneStr+", baseGeneIdx: "+baseGeneIdx);
									NodeList hitHspsList = hitElement.getElementsByTagName("Hit_hsps");
									for (int v = 0; v < hitHspsList.getLength(); v++){
										Node hitHspsNode = hitHspsList.item(v);
										if (hitHspsNode.getNodeType() == Node.ELEMENT_NODE) {
											Element hitHspsElement = (Element) hitHspsNode;
											NodeList hspList = hitHspsElement.getElementsByTagName("Hsp");
				
											for (int w = 0; w < hspList.getLength(); w++) {
												Node hspNode = hspList.item(w);
												if (hspNode.getNodeType() == Node.ELEMENT_NODE) {
													Element hspElement = (Element) hspNode;
													String tmpStr = getTagValue("Hsp_evalue",hspElement);
													double expVal = new BigDecimal(tmpStr).doubleValue();
													if (expVal < minExpVal) {
														minExpVal = expVal;
													}
												}
											}
										}
										if (v == (hitHspsList.getLength()-1)) {
											//System.out.println("tmpCnt: "+tmpCnt);
											MainClass.C[baseOffset+baseGeneIdx][queryOffset+s] = minExpVal;
                                                                                        //tmpCnt++;
										}
									}
								}
									
							}
							//System.out.println("Expectation Value: "+minExpVal);
							if(minExpVal < expValueThreshold)
								MainClass.P[queryOffset + hitsCounter][baseIdx] = 1;
							hitsCounter++;
						}
					}
				}
			} 
			//System.out.println("Number of -filtered- hits: "+hitsCounter);
		} catch (IOException e) {
			System.err.println("Error: " + e.getMessage());
		} catch (javax.xml.parsers.ParserConfigurationException e) {
			System.err.println("Error: " + e.getMessage());
		} catch (org.xml.sax.SAXException e){
			System.err.println("Error: " + e.getMessage());
		}

                boolean success = fXmlFile.delete();
                if(!success)
                    throw new IllegalArgumentException("Delete xml: deletion failed.");
	}
	
	private static String getTagValue(String sTag, Element eElement){
		NodeList nlList= eElement.getElementsByTagName(sTag).item(0).getChildNodes();
		Node nValue = (Node) nlList.item(0); 
		
		return nValue.getNodeValue();    
	}
}
