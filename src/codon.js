// Codon.js

/**
* @namespace Codon
*/
function Codon(){
    // Code goes here.
    
    const DNA_BASES = 'ACGT';
    const RNA_BASES = 'ACGU';
    const quatBases = '0123';
    
    const STOP = 'STOP';
    
    const COMPLEMENTS = {
        A: 'T',
        C: 'G',
        G: 'C',
        T: 'A',
    };
    
    const COMPLEMENTS_RNA = {
        A: 'U',
        C: 'G',
        G: 'C',
        U: 'A',
    };
    
    const AMINO_ACIDS = {
        W: {slc: 'W', tlc: 'Trp', name: 'Tryptophan'},
        C: {slc: 'C', tlc: 'Cys', name: 'Cystine'},
        G: {slc: 'G', tlc: 'Gly', name: 'Glycine'},
        R: {slc: 'R', tlc: 'Arg', name: 'Arginine'},
        S: {slc: 'S', tlc: 'Ser', name: 'Serine'},
        T: {slc: 'T', tlc: 'Thr', name: 'Threonine'},
        A: {slc: 'A', tlc: 'Ala', name: 'Alanine'},
        P: {slc: 'P', tlc: 'Pro', name: 'Proline'},
        F: {slc: 'F', tlc: 'Phe', name: 'Phenylalanine'},
        L: {slc: 'L', tlc: 'Leu', name: 'Leucine'},
        V: {slc: 'V', tlc: 'Val', name: 'Valine'},
        I: {slc: 'I', tlc: 'Ile', name: 'Isoleucine'},
        M: {slc: 'M', tlc: 'Met', name: 'Methionine'},
        Q: {slc: 'Q', tlc: 'Gln', name: 'Glutamine'},
        H: {slc: 'H', tlc: 'His', name: 'Histidine'},
        D: {slc: 'D', tlc: 'Asp', name: 'Aspartic acid'},
        E: {slc: 'E', tlc: 'Glu', name: 'Glutamic acid'},
        K: {slc: 'K', tlc: 'Lys', name: 'Lysine'},
        N: {slc: 'N', tlc: 'Asn', name: 'Asparagine'},
        Y: {slc: 'Y', tlc: 'Tyr', name: 'Tyrosine'},
    };
    
    /**
    * Maps the 64 combinations of codons to their respective amino acids
    */
    const CODONS = {
        UUU: AMINO_ACIDS.F,
        UUC: AMINO_ACIDS.F,
        UUA: AMINO_ACIDS.L,
        UUG: AMINO_ACIDS.L,
        UCU: AMINO_ACIDS.S,
        UCC: AMINO_ACIDS.S,
        UCA: AMINO_ACIDS.S,
        UCG: AMINO_ACIDS.S,
        UAU: AMINO_ACIDS.Y,
        UAC: AMINO_ACIDS.Y,
        UAA: STOP,
        UAG: STOP,
        UGU: AMINO_ACIDS.C,
        UGC: AMINO_ACIDS.C,
        UGA: STOP,
        UGG: AMINO_ACIDS.W,
        CUU: AMINO_ACIDS.L,
        CUC: AMINO_ACIDS.L,
        CUA: AMINO_ACIDS.L,
        CUG: AMINO_ACIDS.L,
        CCU: AMINO_ACIDS.P,
        CCC: AMINO_ACIDS.P,
        CCA: AMINO_ACIDS.P,
        CCG: AMINO_ACIDS.P,
        CAU: AMINO_ACIDS.H,
        CAC: AMINO_ACIDS.H,
        CAA: AMINO_ACIDS.Q,
        CAG: AMINO_ACIDS.Q,
        CGU: AMINO_ACIDS.R,
        CGC: AMINO_ACIDS.R,
        CGA: AMINO_ACIDS.R,
        CGG: AMINO_ACIDS.R,
        AUU: AMINO_ACIDS.I,
        AUC: AMINO_ACIDS.I,
        AUA: AMINO_ACIDS.I,
        AUG: AMINO_ACIDS.M,
        ACU: AMINO_ACIDS.T,
        ACC: AMINO_ACIDS.T,
        ACA: AMINO_ACIDS.T,
        ACG: AMINO_ACIDS.T,
        AAU: AMINO_ACIDS.N,
        AAC: AMINO_ACIDS.N,
        AAA: AMINO_ACIDS.K,
        AAG: AMINO_ACIDS.K,
        AGU: AMINO_ACIDS.S,
        AGC: AMINO_ACIDS.S,
        AGA: AMINO_ACIDS.R,
        AGG: AMINO_ACIDS.R,
        GUU: AMINO_ACIDS.V,
        GUC: AMINO_ACIDS.V,
        GUA: AMINO_ACIDS.V,
        GUG: AMINO_ACIDS.V,
        GCU: AMINO_ACIDS.A,
        GCC: AMINO_ACIDS.A,
        GCA: AMINO_ACIDS.A,
        GCG: AMINO_ACIDS.A,
        GAU: AMINO_ACIDS.D,
        GAC: AMINO_ACIDS.D,
        GAA: AMINO_ACIDS.E,
        GAG: AMINO_ACIDS.E,
        GGU: AMINO_ACIDS.G,
        GGC: AMINO_ACIDS.G,
        GGA: AMINO_ACIDS.G,
        GGG: AMINO_ACIDS.G,
    };
    
    function removeDuplicates(inputArray){
        return Array.from(new Set(inputArray));
    }
    
    /**
    * Generate a quaternary (0-3) numerical sequence from a character sequence
    * @param {string} inputSequence - The character sequence to convert.
    * @param {boolean} rnaFlag - The input sequence is RNA.
    */
    quaternarySequence = function(inputSequence, rnaFlag = false){
        var output = '';
        
        for(var i = 0; i < inputSequence.length; i++){
            output += quatBases.charAt((rnaFlag ? RNA_BASES : DNA_BASES).indexOf(inputSequence[i]));
        }
        
        return output;
    }
    
    /**
    * Generate a random base (ACGT)
    * @param {boolean} rnaFlag - If true, generates Uracil (U) instead of Thymine (T).
    */
    randomBase = function(rnaFlag = false){
        var randIndex = Math.floor(Math.random()*4);
        var base;
        if(rnaFlag){
            base = RNA_BASES.charAt(randIndex);
        }else{
            base = DNA_BASES.charAt(randIndex);
        }
        return base;
    }
    
    /**
    * Count the number of occurences of a base in an input sequence
    * @param {string} inputSequence - The sequence to search
    * @param {string} base - The base to count
    */
    baseCount = function(inputSequence, base){
        var re = new RegExp(base, 'g');
        return inputSequence.match(re || []).length;
    }
    
    
    /**
    * Count the number of occurences of a given k-mer in an input sequence
    * @parm {string} inputSequence - The sequence to search
    * @param {string} kmer - The k-mer to count
    */
    kmerCount = function(inputSequence, kmer){
        var count = 0;
        for(var i = 0; i < inputSequence.length - kmer.length + 1; i++){
            if(inputSequence.substring(i, i+kmer.length) == kmer){                
                count++;
            }
        }
        return count;
    }
    
    /**
    * Find all starting positions in the input sequence where the k-mer appears as a substring
    */
    kmerPositions = function(inputSequence, kmer){
        var positions = [];
        for(var i = 0; i < inputSequence.length - kmer.length + 1; i++){
            if(inputSequence.substring(i, i+kmer.length) == kmer){
                positions.push(i);
            }
        }
        return positions;
    }
    
    /**
    * Find the most frequent k-mers in an input sequence
    * @param {string} inputSequence - the sequence to search
    * @param {number} kmerLength - search for the most frequent k-mers of this length
    * @param {number} minFrequency - only return the most frequent pattern(s) if they occur at least this many times
    */
    kmerMostFrequent = function(inputSequence, kmerLength, minFrequency = 1){
        var mostFrequentPatterns = [];
        var count = [];
        for(var i = 0; i < inputSequence.length - kmerLength + 1; i++){
            var kmer = inputSequence.substring(i, i+kmerLength);
            count[i] = kmerCount(inputSequence, kmer)
        }
        var maxCount = Math.max.apply(null, count);
        
        if(maxCount < minFrequency) {
            return [];
        }
        
        for(var i = 0; i < inputSequence.length - kmerLength + 1; i++){
            if(count[i] == maxCount){
                mostFrequentPatterns.push(inputSequence.substring(i, i+kmerLength));
            }
        }

        return removeDuplicates(mostFrequentPatterns);
    }
    
    /**
    * Find k-mers which occur at a given frequency within a window size
    * @param {string} inputSequence - the sequence to search
    * @param {number} kmerLength - search for k-mers of this length
    * @param {number} frequency - only report k-mers which occur at least this many times within the given window size
    * @param {number} windowSize - the window size to be searched
    */
    findClumps = function(inputSequence, kmerLength, frequency, windowSize){
        var clumps = [];
        for(var i = 0; i < inputSequence.length - windowSize + 1; i++){
            var currentWindow = inputSequence.substring(i, i+windowSize);
            // Merge the output array with the "clumps" array
            Array.prototype.push.apply(clumps, kmerMostFrequent(currentWindow, kmerLength, frequency));
        }
        
        return removeDuplicates(clumps);
    }
    
    /**
    * Generate a random sequence of bases, given an input length
    * @param {number} sequenceLength - The number of bases to generate
    * @param {boolean} rnaFlag - If true, generates a sequence of RNA instead of DNA.
    * @param {boolean} logDetails - If true, prints the details of the sequence to the console after generation.
    */
    randomBaseSequence = function(sequenceLength, rnaFlag = false, logDetails = false){
        
        var sequence = '';       
        
        var numG = 0, numC = 0, numA = 0, numT = 0, numU = 0;
        
        for(var i = 0; i < sequenceLength; i++){
            var randBase = randomBase(rnaFlag);
            sequence += randBase;
            switch(randBase){
                case 'T':
                    numT++;
                    break;
                case 'C':
                    numC++;
                    break;
                case 'A':
                    numA++;
                    break;
                case 'G':
                    numG++;
                    break;
                case 'U':
                    numU++;
                    break;
            }
        }
        
        if(logDetails){
            console.log('Generating random ' + (rnaFlag ? 'RNA' : 'DNA') + ' sequence with '+ sequenceLength +' base pairs: \n');
            
            console.log(sequence);
            
            console.log(
                'A: ' + numA + '\n' +
                'C: ' + numC + '\n' +
                'G: ' + numG + '\n' +
                (rnaFlag ? 'U: ' + numU : 'T: ' + numT) + '\n'
                
            )
        }
        
        return sequence;
    }
    
    /**
    * Generate a single random codon
    * @param {boolean} rnaFlag - If true, generates a codon which may contain Uracil (U) instead of Thymine (T)
    */
    randomCodon = function(rnaFlag = false){
        return randomBaseSequence(3, rnaFlag);
    }
    
    randomCodonSequence = function(numCodons, rnaFlag = false){
        var sequence = '';
        
        for(var i = 0; i < numCodons; i++){
            sequence += randomCodon(rnaFlag);
        }
        
        return sequence;
    }  
    
    /**
    * Calculate the GC content of a given input sequence
    * @param {string} inputSequence - The sequence to analyze
    */
    gcContent = function(inputSequence){
        return ((inputSequence.match(/G/g) || []).length + (inputSequence.match(/C/g) || []).length) / inputSequence.length;
    }
    
    /**
    * Given an input sequence, generate the complementary strand
    */
    complement = function(inputSequence, rnaFlag = false){
        var output = '';
        for(var i = 0; i < inputSequence.length; i++){
            output = (rnaFlag ? COMPLEMENTS_RNA[inputSequence[i]] : COMPLEMENTS[inputSequence[i]]) + output;
        }
        return output;
    }
    
    return{
        'AMINO_ACIDS': AMINO_ACIDS,
        'CODONS': CODONS,
        
        'quaternarySequence': quaternarySequence,
        'baseCount': baseCount,
        'kmerCount': kmerCount,
        'kmerPositions': kmerPositions,
        'kmerMostFrequent': kmerMostFrequent,
        'randomBase': randomBase,
        'randomBaseSequence': randomBaseSequence,
        'randomCodon': randomCodon,
        'randomCodonSequence': randomCodonSequence,
        'gcContent': gcContent,
        //'complement': complement,
    }
    
}