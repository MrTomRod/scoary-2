"use strict"

const colorDict = {
    'g+t+': '#1f78b4',
    'g+t-': '#a6cee3',
    'g+t?': '#e0f4ff',
    'g-t+': '#b2df8a',
    'g-t-': '#33a02c',
    'g-t?': '#edffdf',
}


let _domResolve
const documentReadyPromise = new Promise(function (resolve) {
    _domResolve = resolve
})
document.addEventListener('DOMContentLoaded', _domResolve)


// get trait name
const urlParams = new URLSearchParams(window.location.search)
if (!urlParams.has('trait')) {
    const msg = 'Cannot proceed: this site has to be accessed with an URLSearchParam, i.e. trait.html?trait=<TRAIT>'
    alert(msg)
    throw Error(msg)
}
const trait = urlParams.get('trait')

Papa.execPromise = function (file, config) {
    return new Promise(function (complete, error) {
        config.complete = complete
        config.error = error
        Papa.parse(file, config)
    })
}

// load meta json
const metaPromise = fetch(`traits/${trait}/meta.json`)
    .then(response => response.json())

metaPromise.then(metaData => {
    console.log('metadata-div', metaData)
    document.getElementById('metadata-div').innerText = JSON.stringify(metaData, null, '\t')
})

// Load isolate -> class (-> numeric value)
const valuesPromise = Papa.execPromise(`traits/${trait}/values.tsv`, {
    header: true, download: true, skipEmptyLines: true, delimiter: '\t', newline: '\n',
}).then((valuesData) => {
    valuesData.isNumeric = valuesData.meta.fields.includes('value')
    if (!valuesData.isNumeric) document.getElementById('histogram-container').remove()
    valuesData.isolateToClass = Object.fromEntries(valuesData.data.map(x => [x['isolate'], x['class'].toLowerCase()]))
    if (valuesData.isNumeric) {
        valuesData.isolateToValue = Object.fromEntries(valuesData.data.map(x => [x['isolate'], x['value']]))
    }
    return valuesData
})


// load table
const defaultHiddenCols = ['sensitivity', 'specificity', 'odds_ratio', 'contrasting', 'supporting', 'opposing', 'best', 'worst', 'pval']
const tablePromise = Papa.execPromise(`traits/${trait}/result.tsv`, {
    header: true, download: true, skipEmptyLines: true, delimiter: '\t', newline: '\n',
}).then(tableData => {
    // Preprocess data
    tableData.columns = tableData.meta.fields.map(col => {
        return {'data': col, 'label': col, 'title': col}
    })
    tableData.index = tableData['data'].map(col => col['Gene'])
    tableData.orderCol = tableData.meta.fields.indexOf(
        tableData.meta.fields.includes('pval_empirical') ? 'pval_empirical' : 'qval')
    tableData.hiddenCols = defaultHiddenCols.map(col => tableData.meta.fields.indexOf(col))

    return tableData
})

// load table into DataTable
Promise.all([tablePromise, documentReadyPromise]).then(([tableData, _]) => {
    tableData.table = $('#result-tsv').DataTable({
        data: tableData.data,
        columns: tableData.columns,
        paging: false,
        responsive: true,
        "order": [[tableData.orderCol, "asc"]],
        "columnDefs": [{
            "targets": tableData.hiddenCols, "visible": false, "searchable": false
        }]
    })
    return tableData
}).then((tableData) => {
    // load buttons to toggle columns
    const target = document.getElementById('vis-toggler')
    let htmls = []
    tableData.meta.fields.forEach((col, i) => {
        htmls.push(`<button class="toggle-vis" data-column="${i}">${col}</button>`)
    })
    target.innerHTML = htmls.join(' - ')
    $('.toggle-vis').on('click', function (e) {
        const column = tableData.table.column($(this).attr('data-column'))
        column.visible(!column.visible())
    })
})

// Load coverage matrix data
const coverageMatrixDataPromise = Papa.execPromise(`traits/${trait}/coverage-matrix.tsv`, {
    header: true, download: true, skipEmptyLines: true, delimiter: '\t', newline: '\n',
})

const sanitizeKey = function (key) {
    // remove whitespace and dots
    return encodeURI(key).replace(/\./g, '%2E')
}

const sanitizeData = function (data) {
    data.forEach(function (row) {
        for (const oldKey in row) {
            const newKey = sanitizeKey(oldKey)
            if (oldKey !== newKey) {
                row[newKey] = row[oldKey]
                delete row[oldKey]
            }
        }
    })
    return data
}

// combine values.tsv and coverage_matrix.tsv, sanitize columns
const coverageMatrixPromise = Promise.all([coverageMatrixDataPromise, valuesPromise, metaPromise, documentReadyPromise])
    .then(([cmValues, valuesData, metaData, _]) => {
        // sanitize columns: may not contain dot or whitespace -.-
        const deSanitizeDict = Object.fromEntries(cmValues.meta.fields.map(oldKey => [sanitizeKey(oldKey), oldKey]))
        deSanitizeDict['isolate'] = 'Isolate'
        deSanitizeDict['class'] = 'Class'
        deSanitizeDict['value'] = 'Value'
        cmValues.isNumeric = valuesData.isNumeric
        cmValues.data = sanitizeData(cmValues.data)
        cmValues.meta.fields = cmValues.meta.fields.map(k => sanitizeKey(k))
        const isolates = cmValues.meta.fields.slice(1)  // remove index

        // convert gene list to count matrix
        cmValues.isGeneList = (metaData['genes-content-type'] === 'gene-list')
        if (cmValues.isGeneList) {
            // convert csv to array: 'gene1,gene2' => ['gene1', 'gene2']
            const geneList = cmValues.data
            for (const d of geneList) {
                for (const i of isolates) {
                    d[i] = (d[i] === '') ? [] : d[i].split(',')
                }
            }
            cmValues.geneList = geneList

            // convert gene list to gene  count: ['gene1', 'gene2'] => 2
            const geneCount = []
            for (const d of geneList) {
                const dd = {...d}
                for (const i of isolates) {
                    dd[i] = dd[i].length.toString()
                }
                geneCount.push(dd)
            }

            cmValues.data = geneCount
        }

        // add class column
        cmValues.meta.fields.push('class')
        for (const isolateData of cmValues.data) {
            isolateData['class'] = valuesData.isolateToClass[isolateData['isolate']]
        }

        // add values column
        if (valuesData.isNumeric) {
            cmValues.meta.fields.push('value')
            for (const isolateData of cmValues.data) {
                isolateData['value'] = valuesData.isolateToValue[isolateData['isolate']]
            }
        }

        cmValues.isolates = isolates
        cmValues.deSanitizeDict = deSanitizeDict
        cmValues.index = cmValues.data.map(col => col['isolate'])
        cmValues.columns = cmValues.meta.fields
            .map(col => ({'data': col, 'label': deSanitizeDict[col], 'title': deSanitizeDict[col]}))


        console.log('cmValues', cmValues)
        return cmValues
    })

// plot coverage matrix
const coverageMatrixPlotPromise = coverageMatrixPromise.then((cmValues) => {
    cmValues.table = $('#coverage-matrix').DataTable({
        data: cmValues.data,
        columns: cmValues.columns,
        paging: false,
        searching: false,
        columnDefs: [{
            "targets": [...Array(cmValues.isolates.length + 1).keys()], "orderable": false, // "searchable": false
        }],
        order: [[cmValues.meta.fields.length - 1, "desc"]]  // last col: class or numeric value
    })
    return cmValues
})

const calcLayout = function (layout, heightMultiplier) {
    const plotContainer = document.getElementById('container-right')
    let height = window.innerHeight
    let titleWidth
    for (const titleElement of plotContainer.querySelectorAll('.plot-title')) {
        height -= titleElement.offsetHeight
        titleWidth = titleElement.offsetWidth
    }

    if (window.innerWidth >= 1100) {
        layout.width = titleWidth - 5
        layout.height = height * heightMultiplier - 10
    } else {
        const size = titleWidth * 0.8// square
        layout.width = size
        layout.height = size
    }

    // ensure minimum of 350x350
    layout.width = Math.max(layout.width, 350)
    layout.height = Math.max(layout.height, 350)

    return layout
}

// plot Pie chart for gene
const plotPie = function ({targetId, nameId, heightMultiplier, gene, tableData}) {
    const targetElement = document.getElementById(targetId)
    const targetNameElement = document.getElementById(nameId)
    const geneId = tableData.index.indexOf(gene)
    const colNames = ['g+t+', "g+t-", "g-t+", "g-t-"]

    Plotly.purge(targetElement)
    targetNameElement.innerHTML = gene

    const [gptp, gptn, gntp, gntn] = colNames.map(col => tableData.data[geneId][col])

    const data = [{
        values: [gptp, gptn, gntp, gntn],
        labels: colNames, type: 'pie',
        marker: {colors: colNames.map(col => colorDict[col])}
    }]
    const layout = calcLayout({}, heightMultiplier)
    Plotly.newPlot(targetElement, data, layout)
}

// Plot first gene by default
// Add click listeners to genes
Promise.all([tablePromise, coverageMatrixPlotPromise]).then(([tableData, cmValues]) => {
    const heightMultiplier = cmValues.isNumeric ? 0.5 : 1
    let currentGene
    const plotWrapper = function () {
        plotPie({
            targetId: 'pie',
            nameId: 'pie-gene-name',
            heightMultiplier: heightMultiplier,
            gene: currentGene,
            tableData: tableData
        })
    }

    // set currentGene to first gene in table and draw plot
    currentGene = tableData.index[0]
    plotWrapper()

    $('#coverage-matrix_wrapper thead th:not(:first):not(:last)').click((event) => {
        currentGene = event.target.textContent
        plotWrapper()
    })
    $('#result-tsv tbody tr td:first-child').click((event) => {
        currentGene = event.target.textContent
        plotWrapper()
    })

    let waitUntilResizeOverTimeout
    window.addEventListener('resize', () => {
        clearTimeout(waitUntilResizeOverTimeout);
        waitUntilResizeOverTimeout = setTimeout(plotWrapper, 300);
    })
})

// plot Histogram for gene
const plotHistogram = function ({targetId, nameId, gene, coverageMatrix, valuesData}) {
    const sanitizedGene = sanitizeKey(gene)
    const targetElement = document.getElementById(targetId)
    const targetNameElement = document.getElementById(nameId)

    Plotly.purge(targetElement)
    targetNameElement.innerHTML = gene


    const [gptp, gptn, gptu, gntp, gntn, gntu] = [[], [], [], [], [], []]
    for (const isolateData of coverageMatrix) {
        const isolate = isolateData['isolate']
        const genePositive = isolateData[sanitizedGene].toString() !== '0'
        const traitPositive = valuesData.isolateToClass[isolate]
        const traitValue = valuesData.isolateToValue[isolate]

        if (genePositive) {
            if (traitPositive === 'true') {
                gptp.push(traitValue)
            } else if (traitPositive === 'false') {
                gptn.push(traitValue)
            } else {
                gptu.push(traitValue)
            }
        } else {
            if (traitPositive === 'true') {
                gntp.push(traitValue)
            } else if (traitPositive === 'false') {
                gntn.push(traitValue)
            } else {
                gntu.push(traitValue)
            }
        }
    }

    const data = [
        {name: 'g+t+', x: gptp},
        {name: 'g+t-', x: gptn},
        {name: 'g+t?', x: gptu},
        {name: 'g-t+', x: gntp},
        {name: 'g-t-', x: gntn},
        {name: 'g-t?', x: gntu},
    ]
        .map(trace => {
            trace.type = 'histogram'
            trace.marker = {color: colorDict[trace.name]}
            return trace
        })
        .filter(trace => trace.x.length)  // avoid error messages

    const layout = calcLayout({
        barmode: "stack",
        xaxis: {title: "numerical value"},
        yaxis: {title: "count"}
    }, 0.5)

    Plotly.newPlot(targetElement, data, layout)
}

Promise.all([tablePromise, metaPromise, valuesPromise, coverageMatrixPromise])
    .then(([tableData, metaData, valuesData, coverageMatrixValues]) => {
        if (!valuesData.isNumeric) {
            console.log('no numeric data')
            return
        }

        let currentGene
        const plotWrapper = function () {
            plotHistogram({
                targetId: 'histogram',
                nameId: 'histogram-gene-name',
                gene: tableData.index[0],
                coverageMatrix: coverageMatrixValues.data,
                valuesData: valuesData
            })
        }
        // set currentGene to first gene in table and draw plot
        currentGene = tableData.index[0]
        plotWrapper()

        $('#coverage-matrix_wrapper thead th:not(:first):not(:last)').click((event) => {
            currentGene = event.target.textContent
            plotWrapper()
        })
        $('#result-tsv tbody tr td:first-child').click((event) => {
            currentGene = event.target.textContent
            plotWrapper()
        })

        let waitUntilResizeOverTimeout
        window.addEventListener('resize', () => {
            clearTimeout(waitUntilResizeOverTimeout);
            waitUntilResizeOverTimeout = setTimeout(plotWrapper, 300);
        })
    })
