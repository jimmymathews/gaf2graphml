package com.jm.gaf2graphml;
import uk.ac.ebi.pride.utilities.ols.web.service.client.*;
import uk.ac.ebi.pride.utilities.ols.web.service.model.*;
import uk.ac.ebi.pride.utilities.ols.web.service.config.*;   // import org.springframework.web.client.RestClientException;

import org.apache.commons.io.FileUtils;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.OutputStream;
import java.util.List;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.ArrayList;


import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.OutputKeys;

import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

public class GAF2GraphML 
{
	ArrayList<String> node_categories = new ArrayList<String>();
	ArrayList<String> nodes = new ArrayList<String>();
	ArrayList<String> category_of_corresponding_node = new ArrayList<String>();

	ArrayList<String> rsources = new ArrayList<String>();
	ArrayList<String> rtargets = new ArrayList<String>();
	ArrayList<String> relation_names = new ArrayList<String>();
	ArrayList<String> relation_evidence = new ArrayList<String>();
	ArrayList<String> relation_reference = new ArrayList<String>();

	// boolean evidence_flag  = false;
	// boolean reference_flag = false;
	//* Always write the evidence and reference data into the graphml file, regardless of whether or not gman is ready for them

	String main_filename = "";
	String extension = "";

    OLSWsConfig config;
    OLSClient ols_client;

	Document doc;
	Element gml;
	Element categories_graph;
	Element nodes_graph;

	int color_index = 0;

    public static void main(String[] args)
    {
    	GAF2GraphML g2g = new GAF2GraphML();
    	if(args.length>0)
	    	g2g.convert(args[0]);
	    else
	    	System.out.println("You need to supply a .gaf file as a command-line argument");

	    // g2g.testXMLwriting();
    }

    public GAF2GraphML()
    {
        config = new OLSWsConfig();
        ols_client = new OLSClient(config);
    }

    public void convert(String filename)
    {
		this.setupCategories();
		if(!this.processFilename(filename))
			return;
		this.getDataFromFile(filename);
		this.writeDataToFile(this.main_filename + "." + "graphml");

		System.out.println("Converted GO-flavored gene assocation file "+ filename+ " to GraphML in format readable by gman");
    }

    public boolean processFilename(String filename)
    {
    	Pattern file_pattern = Pattern.compile("(.*)\\.(\\w+)$");
        Matcher m = file_pattern.matcher(filename);
        if (m.find())
        {
            this.main_filename = m.group(1);
            this.extension = m.group(2);
        }
        if(!(this.extension.equals("gaf") || this.extension.equals("GAF")))
        {
        	System.out.println("Not a GAF file");
        	return false;
        }
        return true;
    }

    public void setupCategories()
    {
    	this.addNodeCategory("protein");   // this.add_node_category("evidence");
    	this.addNodeCategory("definition");
    	this.addNodeCategory("process");
    	this.addNodeCategory("function");
    	this.addNodeCategory("component");
    	// this.addNodeCategory("reference");
    	this.addNodeCategory("definition");
    	this.addNodeCategory("other");
    }

    public void addNodeCategory(String cat)
    {
    	for(int i=0; i<node_categories.size(); i++)
    	{
    		if(node_categories.get(i).equals(cat))
    			return;
    	}
    	node_categories.add(cat);
    }

    public void getDataFromFile(String filename)
    {
    	System.out.println("Processing input file:");
    	try
		{
            File f = new File(filename);
            List<String> lines = FileUtils.readLines(f, "UTF-8");
            int count =0;
            int total =0;
            for(String line : lines)
            {
            	total++;
            }
            System.out.println(Integer.toString(total)+ " lines");
            for (String line : lines)
            {
            	count = count + 1;
            	this.updateProgress(count*1.0/(1.0*total), count);
            	this.processLine(line);
	        }
        } catch (IOException e)
        {
            e.printStackTrace();
        }
        System.out.println("Processed each line of "+filename);
    }

    public void processLine(String line)
    {
    	if(line.charAt(0) == '!')
    	{
    		// System.out.println("Ignoring line.");
    		return;
    	}

    	Pattern r = Pattern.compile(this.getLinePattern());
        Matcher m = r.matcher(line);
        if (m.find())
        {
            String protein_name = m.group(3);
            String go_entity_for_involvement = m.group(5);
            String reference_id = m.group(6);
            String evidence_code = m.group(7);
            String aspect = m.group(9);
            String description = m.group(10);
            String object_type = m.group(12);

            // System.out.println("Protein name: "+protein_name);
            // System.out.println("GO entity: "+go_entity_for_involvement);
            // System.out.println("Reference: "+reference_id);
            // System.out.println("Evidence code: "+evidence_code);
            // System.out.println("Description: "+description);

            String relations_string= m.group(16);
			Pattern rels = Pattern.compile("[\\w_]+[(][\\d\\w]{1,30}:[\\w\\d]+[)],?");
        	// Matcher relm = rels.matcher(relations_string);

			List<String> allMatches = new ArrayList<String>();
			Matcher mm = rels.matcher(relations_string);
			while (mm.find())
			{
				allMatches.add(mm.group());
			}
			String[] all_matches = allMatches.toArray(new String[0]);

        	String[] relations = new String[all_matches.length];
        	String[] targets = new String[all_matches.length];

        	for(int j=0; j<all_matches.length; j++)
        	{
        		Pattern rel = Pattern.compile("([\\w_]+)[(]([\\d\\w]{1,30}:[\\w\\d]+)[)]");
        		Matcher relm_single = rel.matcher(all_matches[j]);
        		if(relm_single.find())
        		{
        			String relation = relm_single.group(1);
        			String target = relm_single.group(2);
        			relations[j] = relation;
        			targets[j] = target;
        		}
        	}

        	System.out.println("");
			System.out.println(protein_name + ", "+description + ", \nis associated with "+aspect+"-thing " + go_entity_for_involvement + "\nin the sense of evidence of type "+evidence_code+", according to "+ reference_id);
			if(relations.length>0)
				System.out.println("Also, "+protein_name);
			for(int i=0; i<relations.length; i++)
			{
				System.out.println(relations[i]+ " "+targets[i]);
			}
			System.out.println();


			if(object_type.equals("protein"))
				this.addNode(protein_name, "protein");
			else
			{
				System.out.println("Not a protein?  ("+object_type+")");
				return;
				// this.addNode(protein_name, "other"); //will this happen?
			}

			this.addNode(description, "definition");
			this.addEdge(protein_name, "protein", description, "definition");

			String aspect_string = "";
			if(aspect.equals("P"))
			{
				aspect_string = "process";
			}
			if(aspect.equals("F"))
			{
				aspect_string = "function";
			}
			if(aspect.equals("C"))
			{
				aspect_string = "component";
			}

			// Register the associated entity, and its definition
			String[] lookup = this.lookupEntityNameDefinition(go_entity_for_involvement);
			String name = lookup[0];
			String def = lookup[1];
			this.addNode(name, aspect_string);
			this.addNode(def, "definition");
			this.addEdge(name, aspect_string, def, "definition");

			// Register the association with the associated entity

			ArrayList<String> basic_association_relation_descriptors = new ArrayList<String>();
			basic_association_relation_descriptors.add("associated with");
			basic_association_relation_descriptors.add(evidence_code);
			basic_association_relation_descriptors.add(reference_id);

			this.addEdge(protein_name, "protein", name, aspect_string, basic_association_relation_descriptors);

			// Register the explicit relations to other... proteins? Or GO things, or Ensembl things, or UniProtKB things, or UBERON things, or CL things...
			for(int i=0; i<relations.length; i++)
			{
				lookup = this.lookupEntityNameDefinition(targets[i]);
				name = lookup[0];
				def = lookup[1];
				String type = lookup[2];

				String relation = relations[i];

				ArrayList<String> relation_descriptors = new ArrayList<String>();
				relation_descriptors.add(relation);
				relation_descriptors.add(evidence_code);
				relation_descriptors.add(reference_id);
				if(!type.equals("protein"))
					type = "other";
				this.addNode(name, type);
				this.addNode(def, "definition");
				this.addEdge(name, type, def, "definition");

				this.addEdge(protein_name, "protein", name, type, relation_descriptors);
			}
		}
        else
        {
            System.out.println("NO MATCH: "+line);   //Tested on the basic human database from uniprotkb (october 2017); this exception never occurred (finally)
            System.out.println();
        }
    }

    public String getLinePattern()
    {
		// 1 	DB 	required 	1 	UniProtKB
		// 2 	DB Object ID 	required 	1 	P12345
		// 3 	DB Object Symbol 	required 	1 	PHO3						e.g. protein name
		// 4 	Qualifier 	optional 	0 or greater 	NOT
		// 5 	GO ID 	required 	1 	GO:0003993								GO entity the protein is associated with
		// 6 	DB:Reference (|DB:Reference) 	required 	1 or greater 	PMID:2676709
		// 7 	Evidence Code 	required 	1 	IMP
		// 8 	With (or) From 	optional 	0 or greater 	GO:0000346
		// 9 	Aspect 	required 	1 	F
		// 10 	DB Object Name 	optional 	0 or 1 	Toll-like receptor 4
		// 11 	DB Object Synonym (|Synonym) 	optional 	0 or greater 	hToll|Tollbooth
		// 12 	DB Object Type 	required 	1 	protein
		// 13 	Taxon(|taxon) 	required 	1 or 2 	taxon:9606
		// 14 	Date 	required 	1 	20090118
		// 15 	Assigned By 	required 	1 	SGD
		// 16 	Annotation Extension 	optional 	0 or greater 	part_of(CL:0000576)
		// 17 	Gene Product Form ID 	optional 	0 or 1 	UniProtKB:P12345-2
    	String[] p = new String[18];
    	p[0] = "";
    	p[1] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\<\\: ]+)\\t";
    	p[2] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\<\\: ]+)\\t"; 
    	p[3] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\<\\:/ ]+)\\t";
    	p[4] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\<\\:\\|/ ]*)\\t";
    	p[5] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\<\\:/ ]+)\\t";
    	p[6] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\<:/ ]+)\\t";
    	p[7] = "(\\w{2,3})\\t";
    	p[8] = "([\\w\\d_\\-\\+\\:\\[\\]\\.\\,\\>\\<\\| ]*)\\t";
    	p[9] = "([PFC])\\t";
    	p[10] = "([\\w\\d_\\-\\+\\,\\(\\)\\'\\,\\>\\<\\./\\[\\] \\:]+)\\t";
    	p[11] = "([\\w\\d_\\-\\+\\|/\\[\\]\\.\\,\\>\\<\\'\\* ]*)\\t";
    	p[12] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\< /]+)\\t";
    	p[13] = "([\\w\\d_\\-\\+\\[\\]:\\|\\.\\,\\>\\< ]*)\\t";
    	p[14] = "(\\d{8})\\t";
    	p[15] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\< ]+)\\t";
    	p[16] = "([\\w\\d()_\\-\\+\\[\\]\\:\\|\\,\\>\\<\\. ]*)\\t?";
    	p[17] = "([\\w\\d_\\-\\+\\[\\]\\:\\.\\,\\>\\< ]*)";

       	String pattern = "";
    	for(int i=1; i<18; i++)  //18?
    	{
    		pattern = pattern + p[i];
    	}

    	return pattern;
    }

    public void addNode(String name, String category)
    {
				// 
    	// System.out.println("Adding node  "+name+ "  to category "+category);
    	nodes.add(name);
    	category_of_corresponding_node.add(category);
    }

    public void addEdge(String source, String sourceCat, String target, String targetCat)
    {
    	ArrayList<String> relation_descriptors = new ArrayList<String>();	
    	relation_descriptors.add("");
    	relation_descriptors.add("");
    	relation_descriptors.add("");    	
    	this.addEdge(source, sourceCat, target, targetCat, relation_descriptors);
    }

    public void addEdge(String source, String sourceCat, String target, String targetCat, ArrayList<String> relation_descriptors)
    {
    	// String first_relation = relation_descriptors.get(0);
    	rsources.add(source);   	
    	rtargets.add(target);
    	relation_names.add(relation_descriptors.get(0));
    	relation_evidence.add(relation_descriptors.get(1));
    	relation_reference.add(relation_descriptors.get(2));

    	// System.out.println("Adding edge linking  "+source+ " to "+target + " ("+first_relation+")");

	    //need to decide how to handle annotations of edges... in gman you gave up having categories for edges, but 
    	//that would now mean that all information about an edge must be in its single string description.
    	//so "regulates: " will become "regulates, IDA, PMID:808080: "
    	//that's bad, you need to be able to select them off and on.
    	//*
    	//Solution: Upgrade gman to allow every edge to have a tuple of descriptors. Then 'r'
    	//toggles through them one by one.
    }

	public String[] lookupEntityNameDefinition(String id)
	{
		String id_number = "";
		String database = "";

    	Pattern identifer_pattern = Pattern.compile("(\\w+)\\:([\\w\\d]+)");
        Matcher m = identifer_pattern.matcher(id);
        if (m.find())
        {
            database =  m.group(1);
            id_number = m.group(2);
        }

    	PrintStream original = System.out;
    	System.setOut(new PrintStream(new OutputStream() {
                public void write(int b) {
                    //DO NOTHING
                }
            }));

		String[] l = new String[3];
		String an = "no annotation";
        Term t;
        try
        {	    	
	        t = ols_client.getTermByOBOId(database+":"+id_number,database);
			String def = "";
			String[] de = t.getDescription();
			for(int i=0; i<de.length; i++)
			{
				def = def + ";"+de[i];
			}
			l[0] = t.getName();
			l[1] = def;
			l[2] = "lookedup-type-string";

			if(database.equals("PR") || database.equals("pr"))
			{
				l[2] = "protein";
			}

			// Annotation a = t.getAnnotation();
			// an = a.getValue();
	    }
	    catch(Exception e)
	    {
	    	l[0] = id;
	    	l[1] = "lookup failed "+e.getClass().getCanonicalName();
	    	l[2] = "lookedup-type-string";
	    }
		System.setOut(original);

        // 		throws RestClientException
		// String t.getName();
		// Identifier t.getShortName();
		// Identifier t.getOboId();

		// System.out.println("Label: "+t.getLabel());
		// System.out.println("Name : "+t.getName());

		// System.out.println("Looked up " + id + " ... " + l[0]);   /////

		// System.out.println("      Annotation: "+an);
		return l;
        // Term t = cl.getTermByOBOId("GO:0005515","GO");
			        // throws RestClientException

		// ...
		// String t.getName()
		// Identifier t.getShortName()
		// Identifier t.getOboId()

		// System.out.println("Label: "+t.getLabel());
		// System.out.println("Name : "+t.getName());
		// String[] de = t.getDescription();
		// for(int i=0; i<de.length; i++)
		// {
		// 	System.out.println(de[i]);
		// }

	}

	public void writeDataToFile(String filename)
	{
		this.startGraphMLDocument();
		this.addCategories();
		this.addNodes();
		this.addLinks();
		this.writeOut(filename);
	}

	public void startGraphMLDocument()
	{
		try{
			DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
			doc = docBuilder.newDocument();

			gml = doc.createElement("graphml");
			doc.appendChild(gml);

			categories_graph = doc.createElement("graph");
			categories_graph.setAttribute("id","C");
			categories_graph.setAttribute("edgedefault","directed");
			gml.appendChild(categories_graph);

			nodes_graph = doc.createElement("graph");
			nodes_graph.setAttribute("id","N");
			nodes_graph.setAttribute("edgedefault","directed");
			gml.appendChild(nodes_graph);
		}
		catch(ParserConfigurationException pce)
		{
			pce.printStackTrace();
		}
	}

	public void addCategories()
	{
		for(int i=0; i<node_categories.size(); i++)
		{
			this.addCategory(node_categories.get(i), "c"+Integer.toString(i));
		}
	}

	public void addCategory(String cat, String id)
	{
		Element xc = doc.createElement("node");
		xc.setAttribute("id", id);
		categories_graph.appendChild(xc);

		Element d = doc.createElement("data");
		d.setAttribute("key","d1");
		d.appendChild(doc.createTextNode(cat));
		xc.appendChild(d);

		Element cl = doc.createElement("data");
		cl.setAttribute("key","d2");
		cl.appendChild(doc.createTextNode(Integer.toString(this.getNextColorIndex())));
		xc.appendChild(cl);
	}

	public int getNextColorIndex()
	{
		color_index = color_index + 1;
		if(color_index > 7)
			color_index = 0;
		return color_index;
	}

	public void addNodes()
	{
		for(int i=0; i<nodes.size(); i++)
		{
			Element n = doc.createElement("node");
			n.setAttribute("id","n"+Integer.toString(i));
			nodes_graph.appendChild(n);

			Element d = doc.createElement("data");
			d.setAttribute("key","d0");
			d.appendChild(doc.createTextNode(nodes.get(i)));
			n.appendChild(d);

			Element cat = doc.createElement("data");
			cat.setAttribute("key","d1");
			cat.appendChild(doc.createTextNode(category_of_corresponding_node.get(i)));
			n.appendChild(cat);
		}
	}

	public void addLinks()
	{
	// ArrayList<String> rsources = new ArrayList<String>();
	// ArrayList<String> rtargets = new ArrayList<String>();
	// ArrayList<String> relation_names = new ArrayList<String>();
	// ArrayList<String> relation_evidence = new ArrayList<String>();
	// ArrayList<String> relation_reference = new ArrayList<String>();

	// 	pugi::xml_node xn = nodes_graph.append_child("edge");
	// 	xn.append_attribute("id") = (n->get_id() + l->get_end_node()->get_id()).c_str();
	// 	xn.append_attribute("source") = n->get_id().c_str();
	// 	xn.append_attribute("target") = l->get_end_node()->get_id().c_str();

	// 	pugi::xml_node c = xn.append_child("data");
	// 	c.append_attribute("key") = "d3";
	// 	c.append_child(pugi::node_pcdata).set_value(l->get_name().c_str());				

		for(int i=0; i<relation_names.size(); i++)
		{
			String source_id = this.lookupNodeID(rsources.get(i));
			String target_id = this.lookupNodeID(rtargets.get(i));
			Element e = doc.createElement("edge");
			e.setAttribute("id",source_id + target_id);
			e.setAttribute("source",source_id);
			e.setAttribute("target",target_id);
			nodes_graph.appendChild(e);

			Element r1 = doc.createElement("data");
			r1.setAttribute("key","d3");
			r1.appendChild(doc.createTextNode(relation_names.get(i)));
			e.appendChild(r1);			

			Element r2 = doc.createElement("data");
			r2.setAttribute("key","d4");
			r2.appendChild(doc.createTextNode(relation_evidence.get(i)));
			e.appendChild(r2);			

			Element r3 = doc.createElement("data");
			r3.setAttribute("key","d5");
			r3.appendChild(doc.createTextNode(relation_reference.get(i)));
			e.appendChild(r3);			
		}
	}

	String lookupNodeID(String contents)
	{
		for(int i=0; i<nodes.size(); i++)
		{
			if(nodes.get(i).equals(contents))
				return "n"+Integer.toString(i);
		}
		System.out.println("  -- Node not found: " + contents );
		return "n-1";
	}

	public void writeOut(String filename)
	{
		try{
			// write the content into xml file
			TransformerFactory transformerFactory = TransformerFactory.newInstance();
			Transformer transformer = transformerFactory.newTransformer();
			DOMSource source = new DOMSource(doc);
			StreamResult result = new StreamResult(new File(filename));

			// Output to console for testing
			// StreamResult result = new StreamResult(System.out);

			transformer.setOutputProperty(OutputKeys.INDENT, "yes");
			transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
			transformer.transform(source, result);

			// System.out.println("File saved.");
		}
		catch(TransformerException tfe)
		{
			tfe.printStackTrace();
		}
	}

    public void testXMLwriting()
    {
		try
		{
			DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder docBuilder = docFactory.newDocumentBuilder();

			// root elements
			Document doc = docBuilder.newDocument();
			Element rootElement = doc.createElement("company");
			doc.appendChild(rootElement);

			// staff elements
			Element staff = doc.createElement("Staff");
			rootElement.appendChild(staff);

			// set attribute to staff element
			Attr attr = doc.createAttribute("id");
			attr.setValue("1");
			staff.setAttributeNode(attr);

			// shorten way
			// staff.setAttribute("id", "1");

			// firstname elements
			Element firstname = doc.createElement("firstname");
			firstname.appendChild(doc.createTextNode("yong"));
			staff.appendChild(firstname);

			// lastname elements
			Element lastname = doc.createElement("lastname");
			lastname.appendChild(doc.createTextNode("mook kim"));
			staff.appendChild(lastname);

			// nickname elements
			Element nickname = doc.createElement("nickname");
			nickname.appendChild(doc.createTextNode("mkyong"));
			staff.appendChild(nickname);

			// salary elements
			Element salary = doc.createElement("salary");
			salary.appendChild(doc.createTextNode("100000"));
			staff.appendChild(salary);

			// write the content into xml file
			TransformerFactory transformerFactory = TransformerFactory.newInstance();
			Transformer transformer = transformerFactory.newTransformer();
			DOMSource source = new DOMSource(doc);
			StreamResult result = new StreamResult(new File("file.xml"));

			// Output to console for testing
			// StreamResult result = new StreamResult(System.out);

			transformer.transform(source, result);

			// System.out.println("File saved!");

		} catch (ParserConfigurationException pce)
		{
			pce.printStackTrace();
		} catch (TransformerException tfe)
		{
			tfe.printStackTrace();
		}
    }


	void updateProgress(double progressPercentage, int count)
	{
		int width = 90; // progress bar width in chars

		// System.out.print("\r[");
		System.out.print("[");
		int i = 0;
		for (; i <= (int)(progressPercentage*width); i++)
		{
			System.out.print(".");
		}
		for (; i < width; i++) {
		  System.out.print(" ");
		}
		System.out.print("] "+Integer.toString(count));
		System.out.println("");
	}
}

        // System.out.println(t.getIri());
		// System.out.println();
		// String[] de = t.getDescription();
		// for(int i=0; i<de.length; i++)
		// {
		// 	System.out.println(de[i]);
		// }
		// System.out.println();

		// System.out.println(t.getOntologyName());
		// System.out.println(t.getScore());
		// System.out.println(t.getOntologyPrefix());
		// System.out.println(t.getOntologyIri());

		// Identifier t.getShortForm()
		// Identifier t.getTermOBOId()
		// Annotation t.getAnnotation()
		// Identifier t.getGlobalId()




			// Element contents_key = doc.createElement("key");
			// contents_key.setAttribute("id","d0");
			// contents_key.setAttribute("for","node");
			// contents_key.setAttribute("attr.name","contents");
			// contents_key.setAttribute("attr.type","string");	
			// doc.appendChild(contents_key);

			// Element category_key = doc.createElement("key");
			// category_key.setAttribute("id","d1");
			// category_key.setAttribute("for","node");
			// category_key.setAttribute("attr.name","category");
			// category_key.setAttribute("attr.type","string");	
			// doc.appendChild(category_key);

			// Element color_key = doc.createElement("key");
			// color_key.setAttribute("id","d2");
			// color_key.setAttribute("for","node");
			// color_key.setAttribute("attr.name","color");
			// color_key.setAttribute("attr.type","int");	
			// doc.appendChild(color_key);

			// Element relation_key = doc.createElement("key");
			// relation_key.setAttribute("id","d3");
			// relation_key.setAttribute("for","edge");
			// relation_key.setAttribute("attr.name","relation");
			// relation_key.setAttribute("attr.type","string");	
			// doc.appendChild(relation_key);
