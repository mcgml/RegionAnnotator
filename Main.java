package nhs.genetics.cardiff;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());

    public static void main(String[] args) {

        if (args.length != 3) {
            System.err.println("Usage: <BedFile> <VepFile> <Transcripts>");
            System.exit(1);
        }

        ArrayList<String> nameTemp = new ArrayList();
        File bedFilePath = new File(args[0]);
        File vepFilePath = new File(args[1]);
        File transcriptFilePath = new File(args[2]);

        String[] filename = bedFilePath.getName().split("\\.");

        //read preferred transcripts
        HashSet<String> preferredTranscripts = new HashSet<>();
        try (BufferedReader transcriptReader = new BufferedReader(new FileReader(transcriptFilePath))) {

            String line;

            while ((line = transcriptReader.readLine()) != null) {
                if (!line.equals("")){
                    String[] fields = line.split("\t");
                    preferredTranscripts.add(fields[1]);
                }
            }

        } catch (IOException e) {
            log.log(Level.SEVERE, e.getMessage());
        }

        //read annotation file
        VEPAnnotationFile annotationFile = new VEPAnnotationFile(vepFilePath);
        annotationFile.parseVEP();

        //read bed file for annotation
        try (PrintWriter printWriter = new PrintWriter(new File(filename[0] + "_Annotated.bed"))){

            try (AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(bedFilePath.toString(), new BEDCodec(BEDCodec.StartOffset.ZERO), false)){
                Iterable<BEDFeature> iter = reader.iterator();

                //loop over bed records
                for (BEDFeature feature : iter) {

                    printWriter.write(feature.getContig() + "\t");
                    printWriter.write(feature.getStart() + "\t");
                    printWriter.write(feature.getEnd() + "\t");

                    if (!feature.getName().equals("")) {
                        printWriter.write(feature.getName());
                        printWriter.write(";");
                    }

                    //loop over annotations for this region & print HGVS
                    if (annotationFile.getTranscriptLevelRecords().containsKey(feature.getContig() + ":" + feature.getStart() + "-" + feature.getEnd())){ //check if annotations exist

                        //loop over annotations
                        for (VEPTranscriptAnnotation ann : annotationFile.getTranscriptLevelRecords().get(feature.getContig() + ":" + feature.getStart() + "-" + feature.getEnd()))
                        {
                            if (preferredTranscripts.size() != 0){ //if transcripts are provided; filter

                                if (preferredTranscripts.contains(ann.getFeature())){
                                    String[] fields = ann.getHgvsCoding().split("del");
                                    nameTemp.add(ann.getFeature() + ":" + ann.getSymbol() + ":" + fields[0]);
                                }

                            } else {
                                String[] fields = ann.getHgvsCoding().split("del");
                                nameTemp.add(ann.getFeature() + ":" + ann.getSymbol() + ":" + fields[0]);
                            }

                        }

                    }

                    //print name field
                    for (int n = 0; n < nameTemp.size(); ++n){
                        if (n != nameTemp.size() - 1){
                            printWriter.write(nameTemp.get(n) + ";");
                        } else {
                            printWriter.write(nameTemp.get(n));
                        }
                    }

                    printWriter.println();
                    nameTemp.clear();
                }

                reader.close();
            } catch (IOException e){
                log.log(Level.SEVERE, e.toString());
            }

            printWriter.close();
        } catch (IOException e){
            log.log(Level.SEVERE, e.toString());
        }

    }
}
