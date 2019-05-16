var codon = new Codon();

function colorBases() {
    $('.bases').each(function () {
        var formatted = $(this).html();
        formatted = formatted.replace(/A/g, '<span class="red">A</span>');
        formatted = formatted.replace(/C/g, '<span class="green">C</span>');
        formatted = formatted.replace(/G/g, '<span class="blue">G</span>');
        formatted = formatted.replace(/T/g, '<span class="orange">T</span>');
        formatted = formatted.replace(/U/g, '<span class="orange">U</span>');
        $(this).html(formatted);
    });
}

function makeRandomSequence() {
    var generatedSequence = randomBaseSequence(33, true, true);
    $('#random-seq').html(generatedSequence);

    updateData(generatedSequence);
}

function handleFileUpload() {
    var file = document.getElementById('file-upload').files[0];

    var filename = file.name;
    var reader = new FileReader();

    reader.onload = function (theFile) {
        console.log('File read of ' + filename + ' complete.');
        var rawText = theFile.target.result;
        var genomeString = parseFasta(rawText, 0, 50);
        $('#loaded-seq').html(genomeString);
        $('#loaded-filename').html(filename);

        updateData(genomeString);
    }

    reader.readAsText(file);
}

function loadFile(filePath) {
    var result = null;
    var xmlhttp = new XMLHttpRequest();
    xmlhttp.open("GET", filePath, false);
    xmlhttp.send();
    if (xmlhttp.status == 200) {
        result = xmlhttp.responseText;
    }
    return result;
}

function loadInitSequence() {
    var file = loadFile('genomes/mycoplasma_genitalium.fa');
    var genomeString = parseFasta(file, 0, 50);

    $('#loaded-seq').html(genomeString);
}

// Update the various data displays based on an input sequence.
function updateData(inputSequence) {
    var rnaFlag = isRNA(inputSequence);

    var mostFrequent = kmerMostFrequent(inputSequence, 4);
    $('#most-frequent').html(mostFrequent.toString());
    $('#most-frequent-copy').html(mostFrequent[0]);

    $('#gc-content').html(Math.floor(gcContent(inputSequence) * 100) + "%");
    $('#complement').html(complement(inputSequence, rnaFlag));

    var acids = '';
    var acidSequence = inputSequence;
    if (!rnaFlag) {
        acidSequence = dnaToRna(inputSequence);
    }
    console.log(acidSequence);
    for (var i = 0; i < acidSequence.length - 2; i += 3) {
        var current = acidSequence.substring(i, i + 3);
        var name = codon.CODONS[current];
        if (name != 'STOP') {
            acids += name.tlc;
        } else {
            acids += name;
        }
    }

    $('#acids').html(acids);

    $('#positions').html(kmerPositions(inputSequence, mostFrequent[0]).toString());

    $('#skews').html(minimumSkew(inputSequence).toString());

    colorBases();
}

loadInitSequence();
makeRandomSequence();
colorBases();