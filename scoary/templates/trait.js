"use strict"


/**
 * String formatting function (similar to Python 'str.format')
 *
 * Example:
 *  - Input: "My Name is {name} and I'm {age} years old!".format({"name": "Jim", "age": "28"})
 *  - Output: "My Name is Jim and I'm 28 years old!"
 */
String.prototype.format = function (dictionary) {
    return this.replace(/\{[0-9a-zA-Z_-]+\}/g, function (prekey) {
        const key = prekey.slice(1, -1)  // remove {}
        const val = dictionary[key]
        if (typeof val !== 'string') {
            throw {
                name: 'StringFormattingException',
                message: `Failed to format string. key=${key} not in dictionary=${JSON.stringify(dictionary)}`,
                toString: function () {
                    return `${this.name}: ${this.message}`
                }
            }
        }
        return val
    })
}


/**
 * Load the name of the trait from URL
 */
const urlParams = new URLSearchParams(window.location.search)
if (!urlParams.has('trait')) {
    const msg = 'Cannot proceed: this site has to be accessed with an URLSearchParam, i.e. trait.html?trait=<TRAIT>'
    alert(msg)
    throw Error(msg)
}
const trait = urlParams.get('trait')


/**
 * Returns a promise for Papa.parse.
 *
 * @param {string} file File to be loaded using Papa.parse.
 * @param {Object} config Configuration object, see https://www.papaparse.com/docs#config.
 */
Papa.execPromise = function (file, config) {
    return new Promise(function (complete, error) {
        config.complete = complete
        config.error = error
        Papa.parse(file, config)
    })
}


/**
 * Remove whitespace and characters like - or . from the key.
 *
 * Required because DataTables otherwise throws an error (Requested unknown parameter)
 * https://datatables.net/manual/tech-notes/4
 *
 * @param {string} key String to be sanitized.
 */
const sanitizeKey = function (key) {
    return encodeURI(key).replace(/\./g, '%2E')
}


/**
 * Sanitize the column names of a table.
 *
 * @param {Array} data Array of columns: [{Isolate: str, orthogene1: any}].
 */
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


/**
 * Create a promise for when the document is ready.
 * Will consistently use promises to manage async operations.
 */
let _domResolve
const documentReadyPromise = new Promise(function (resolve) {
    _domResolve = resolve
})
document.addEventListener('DOMContentLoaded', _domResolve)


/**
 * Promise for config.json. Returns the dictionary object.
 */
const configPromise = fetch(`config.json`)
    .then(response => response.json())
    .then(config => {
        if (config['sanitize-genes']) {
            config.geneListToArray = (d, i) => d[i].split(',').map((gene) => gene.split('|').slice(-1)[0])
        } else {
            config.geneListToArray = (d, i) => d[i].split(',')
        }
        return config
    })


/**
 * Promise for meta.json. Returns the dictionary object.
 */
const metaPromise = fetch(`traits/${trait}/meta.json`)
    .then(response => response.json())


/**
 * Show metadata.
 */
metaPromise.then(metaData => {
    document.getElementById('metadata-div').innerText = JSON.stringify(metaData, null, '\t')
})


/**
 * Promise for values.tsv. Returns the valuesData object.
 *
 * valuesData contains:
 *  - valeusData.data                [{isolate: string, class: "True"|"False"|"", value: float}]
 *  - valeusData.meta                {fields: [str]}
 *  - valuesData.isNumeric           if false, the data is boolean
 *  - valuesData.isolateToTrait      {"isolate": boolean}
 *  - valuesData.isolateToValue      {"isolate": float}  (only if isNumeric)
 *  - valuesData.isolateToNormValue  {"isolate": float(0-1)}  (only if isNumeric)
 *  - valuesData.colorScale          chroma color scale function: colorScale(float) -> "color"  (only if isNumeric)
 */
const valuesPromise = configPromise.then(config => {
    return Papa.execPromise(`traits/${trait}/values.tsv`, {
        header: true, download: true, skipEmptyLines: true, delimiter: '\t', newline: '\n',
    }).then((valuesData) => {
        valuesData.isNumeric = valuesData.meta.fields.includes('value')
        if (!valuesData.isNumeric) document.getElementById('histogram-container').remove()
        valuesData.isolateToTrait = Object.fromEntries(valuesData.data.map(x => [x['isolate'], x['class'].toLowerCase()]))
        if (valuesData.isNumeric) {
            valuesData.isolateToValue = Object.fromEntries(valuesData.data.map(x => [x['isolate'], x['value']]))
            const vals = valuesData.data.map(x => x['value'])
            const valMax = Math.max(...vals)
            const valMin = Math.min(...vals)
            const valDiff = valMax - valMin
            valuesData.isolateToNormValue = Object.fromEntries(valuesData.data.map(
                x => [x['isolate'], (x['value'] - valMin) / valDiff]
            ))
            valuesData.colorScale = chroma.scale(config['color-scale'])
        }

        return valuesData
    })

})


/**
 * Promise for values.tsv. Returns the valuesData object.
 *
 * valuesData contains:
 *  - valuesData.data        [{Gene: string, g+t+: int, ...}]
 *  - valuesData.meta        {fields: [str]}
 *  - valuesData.index       [str]  (list of genes)
 *  - valuesData.columns     [{data: str, label: str, title:str}]  (DataTables format)
 *  - valuesData.hiddenCols  [int]  (indices of hidden columns)
 *  - valuesData.orderCol    int  (index of order column)
 *  - valuesData.table       DataTables object  (is added later, see genesTablePromise)
 */
const tablePromise = configPromise.then(config => {
    return Papa.execPromise(`traits/${trait}/result.tsv`, {
        header: true, download: true, skipEmptyLines: true, delimiter: '\t', newline: '\n',
    }).then(tableData => {
        // Preprocess data
        tableData.columns = tableData.meta.fields.map(col => {
            return {'data': col, 'label': col, 'title': col}
        })
        tableData.index = tableData['data'].map(col => col['Gene'])
        tableData.orderCol = tableData.meta.fields.indexOf(
            tableData.meta.fields.includes('pval_empirical') ? 'pval_empirical' : 'qval')
        tableData.hiddenCols = config['default-hidden-cols'].map(col => tableData.meta.fields.indexOf(col))

        return tableData
    })
})


/**
 * Load Genes table.
 */
const genesTablePromise = Promise.all([tablePromise, documentReadyPromise]).then(([tableData, _]) => {
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


/**
 * Promise for isolate_info.tsv. Returns the isolateInfo object if the file exists.
 *
 * isolateInfo contains:
 *  - isolateInfo.data        [{Isolate: string, ...}]
 *  - isolateInfo.meta        {fields: [str]}
 *  - isolateInfo.index       [str]  (list of genes)
 *  - isolateInfo.columns     [{Isolate: str, ...}]  (DataTables format)
 */
const isolateInfoPromise = Papa.execPromise('isolate_info.tsv', {
    header: true, download: true, skipEmptyLines: true, delimiter: '\t', newline: '\n',
}).then(isolateInfo => {
    // Preprocess data
    isolateInfo.columns = isolateInfo.meta.fields.map(col => {
        return {'data': col, 'label': col, 'title': col}
    })
    isolateInfo.index = isolateInfo['data'].map(col => col['Isolate'])

    return isolateInfo
}).catch(() => {
    console.log('no isolate info')
    return 'no isolate info'
})


/**
 * Promise for coverage-matrix.tsv. Returns the raw output of Papa.parse.
 */
const coverageMatrixDataPromise = Papa.execPromise(`traits/${trait}/coverage-matrix.tsv`, {
    header: true, download: true, skipEmptyLines: true, delimiter: '\t', newline: '\n',
})


/**
 * Promise for coverage-matrix.tsv. Returns the cmValues object.
 *
 * cmValues contains:
 *  - cmValues.data            [{Isolate: string, orthogene1: float, gene2:float, ...}]  (gene count table)
 *  - cmValues.meta            {fields: [str]}
 *  - cmValues.index           [str]  (list of isolates)
 *  - cmValues.isNumeric       if false, gene count data only
 *  - cmValues.deSanitizeDict  {str: str}  (sanitized -> original column names)
 *  - cmValues.sanitizedGenes  [str]  (sanitized orthogene names)
 *  - cmValues.rawGenes        [str]  (raw orthogene names)
 *  - cmValues.isGeneList      boolean
 *  - cmValues.geneList        [{Isolate: str, orthogene1: [genes], ...}]  (list of genes table)
 *  - cmValues.table           DataTables object  (is added later, see coverageMatrixTablePromise)
 */
const coverageMatrixPromise = Promise.all(
    [coverageMatrixDataPromise, valuesPromise, metaPromise, isolateInfoPromise, configPromise, documentReadyPromise]
).then(([cmValues, valuesData, metaData, isolateInfo, config, _]) => {
    // sanitize columns: may not contain dot or whitespace -.-
    const deSanitizeDict = Object.fromEntries(cmValues.meta.fields.map(oldKey => [sanitizeKey(oldKey), oldKey]))
    deSanitizeDict['Isolate'] = 'Isolate'
    deSanitizeDict['class'] = 'Class'
    deSanitizeDict['value'] = 'Value'
    cmValues.isNumeric = valuesData.isNumeric
    cmValues.data = sanitizeData(cmValues.data)
    const rawGenes = cmValues.meta.fields.slice(1)
    cmValues.meta.fields = cmValues.meta.fields.map(k => sanitizeKey(k))
    const sanitizedGenes = cmValues.meta.fields.slice(1)

    // convert gene list to count matrix
    cmValues.isGeneList = (metaData['genes-content-type'] === 'gene-list')
    if (cmValues.isGeneList) {
        // convert csv to array: 'orthogene1,orthogene2' => ['orthogene1', 'orthogene2']
        const geneList = cmValues.data

        for (const d of geneList) {
            for (const i of sanitizedGenes) {
                d[i] = (d[i] === '') ? [] : config.geneListToArray(d, i)
            }
        }

        cmValues.geneList = geneList

        // convert gene list to gene  count: ['orthogene1', 'orthogene2'] => 2
        const geneCount = []
        for (const d of geneList) {
            const dd = {...d}
            for (const i of sanitizedGenes) {
                dd[i] = dd[i].length.toString()
            }
            geneCount.push(dd)
        }

        cmValues.data = geneCount
    }

    // add isolate information
    if (isolateInfo !== 'no isolate info') {
        for (const col of isolateInfo.meta.fields.slice(1)) {
            const sanitizedCol = sanitizeKey(col)
            deSanitizeDict[sanitizedCol] = col
            cmValues.meta.fields.push(sanitizedCol)
            for (const isolateData of cmValues.data) {
                const isolateIndex = isolateInfo.index.indexOf(isolateData['Isolate'])
                if (isolateIndex === -1) {
                    isolateData[sanitizedCol] = ''// no data for this isolate in the table
                } else {
                    isolateData[sanitizedCol] = isolateInfo.data[isolateIndex][col]
                }
            }
        }
    }

    // add class column
    cmValues.meta.fields.push('class')
    for (const isolateData of cmValues.data) {
        isolateData['class'] = valuesData.isolateToTrait[isolateData['Isolate']]
    }

    // add values column
    if (valuesData.isNumeric) {
        cmValues.meta.fields.push('value')
        for (const isolateData of cmValues.data) {
            isolateData['value'] = valuesData.isolateToValue[isolateData['Isolate']]
        }
    }

    // save some variables to promise output
    cmValues.rawGenes = rawGenes
    cmValues.sanitizedGenes = sanitizedGenes
    cmValues.deSanitizeDict = deSanitizeDict
    cmValues.index = cmValues.data.map(col => col['Isolate'])

    return cmValues
})


/**
 * Hide all shown popovers if a click on a non-popover is registered.
 */
let popovers = []
$(document).on("click", function (event) {
    // ignore clicks on .popover or .has-popover
    if (Boolean(event.target.closest('.popover, .has-popover'))) {
        return
    }

    // hide all popovers
    while (popovers.length > 0) {
        const popover = popovers.pop()
        bootstrap.Popover.getInstance(popover).hide()
    }

})


/**
 * Load the coverage matrix table, make it interactive.
 */

const coverageMatrixTablePromise = Promise.all(
    [coverageMatrixPromise, configPromise]
).then(([cmValues, config]) => {
    // Load DataTable
    const generateColumnDef = function (col, i) {
        const colDef = {
            'data': col,
            'label': cmValues.deSanitizeDict[col],
            'title': cmValues.deSanitizeDict[col]
        }
        if (i > 0 && i <= cmValues.sanitizedGenes.length) {
            colDef['className'] = 'gene-count'
        }
        return colDef

    }
    cmValues.table = $('#coverage-matrix').DataTable({
        data: cmValues.data,
        columns: cmValues.meta.fields
            .map((col, i) => generateColumnDef(col, i)),
        paging: false,
        searching: false,
        columnDefs: [{
            "targets": [...Array(cmValues.sanitizedGenes.length + 1).keys()], "orderable": false, // "searchable": false
        }],
        order: [[cmValues.meta.fields.length - 1, "desc"]]  // last col: class or numeric value
    })


    /**
     * Generate title of popover.
     */
    const genePopoverTitle = function () {
        this.classList.add('has-popover')
        popovers.push(this)
        const isolate = this.parentElement.firstChild.textContent
        const geneCoord = this.cellIndex - 1
        const geneSanitized = cmValues.sanitizedGenes[geneCoord]
        const gene = cmValues.deSanitizeDict[geneSanitized]
        return `${gene} x ${isolate}`
    }


    /**
     * Make genes clickable depending on config.json.
     */
    const createGeneElement = function (gene) {
        let geneLink = config['gene-link']
        if (geneLink) {
            geneLink = String(geneLink).format({gene: gene})
            return `<a href='${geneLink}'><code>${gene}</code></a>`
        } else {
            return `<code>${gene}</code>`
        }
    }


    /**
     * Generate content of popover.
     */
    const genePopoverContent = function () {
        const isolate = this.parentElement.firstChild.textContent
        const geneCoord = this.cellIndex - 1
        const geneSanitized = cmValues.sanitizedGenes[geneCoord]
        const orthoGene = cmValues.deSanitizeDict[geneSanitized]
        const isolateId = cmValues.index.indexOf(isolate)
        const isolateData = cmValues.geneList[isolateId]
        let genes = isolateData[geneSanitized]

        if (genes.length) {
            const htmls = genes.map((gene) => createGeneElement(gene))
            return htmls.join('<br>')
        } else {
            return 'no genes'
        }
    }

    // Create popover events using bootstrap (https://getbootstrap.com/docs/5.1/components/popovers/)
    const allGeneCells = document.querySelectorAll('#coverage-matrix td.gene-count')
    Array.from(allGeneCells).forEach((patch, index) => {
        new bootstrap.Popover(patch, {
            container: 'body',
            trigger: 'click',
            html: true,
            title: genePopoverTitle,
            content: genePopoverContent
        })
    })

    return cmValues
})


/**
 * Calculate the width and height of pie plot / histogram based on size of #container-right.
 *
 * Adds these parameters to layout:
 *  - layout.height: int
 *  - layout.width: int
 *
 * @param {Object} layout Configuration for Plotly plot.
 * @param {float|int} heightMultiplier The height will be scaled by this value
 */
const calcLayout = function (layout, heightMultiplier) {
    const plotContainer = document.getElementById('container-right')
    let height = window.innerHeight
    let titleWidth
    for (const titleElement of plotContainer.querySelectorAll('.plot-title')) {
        height -= titleElement.offsetHeight
        titleWidth = titleElement.offsetWidth
    }

    if (window.innerWidth >= 1100) {
        // show plots on the right side, see also trait.css
        layout.width = titleWidth - 5
        layout.height = height * heightMultiplier - 10
    } else {
        // show plots below
        const size = titleWidth * 0.8 // square
        layout.width = size
        layout.height = size
    }

    // ensure minimum of 350x350
    layout.width = Math.max(layout.width, 350)
    layout.height = Math.max(layout.height, 350)

    return layout
}


/**
 * For a gene, calculate the numeric values per category. (Raw data for histogram.)
 *
 * @param {string} gene Name of the orthgene.
 * @param {[]} coverageMatrix Gene count table
 * @param {[]} valuesData Boolean trait table
 */
const calcHistValues = function (gene, coverageMatrix, valuesData) {
    const sanitizedGene = sanitizeKey(gene)
    const [gptp, gptn, gptu, gntp, gntn, gntu] = [[], [], [], [], [], []]
    for (const isolateData of coverageMatrix) {
        const isolate = isolateData['Isolate']
        const genePositive = isolateData[sanitizedGene].toString() !== '0'
        const traitPositive = valuesData.isolateToTrait[isolate]
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
    return [gptp, gptn, gptu, gntp, gntn, gntu]
}


/**
 * For a gene, calculate the number of isolates per category. (Raw data for histogram.)
 *
 * @param {string} gene Name of the orthgene.
 * @param {[]} coverageMatrix Gene count table
 * @param {[]} valuesData Boolean trait table
 */
const calcPieValues = function (gene, coverageMatrix, valuesData) {
    const sanitizedGene = sanitizeKey(gene)
    let [gptp, gptn, gptu, gntp, gntn, gntu] = [0, 0, 0, 0, 0, 0]
    for (const isolateData of coverageMatrix) {
        const isolate = isolateData['Isolate']
        const genePositive = isolateData[sanitizedGene].toString() !== '0'
        const traitPositive = valuesData.isolateToTrait[isolate]

        if (genePositive) {
            if (traitPositive === 'true') {
                gptp += 1
            } else if (traitPositive === 'false') {
                gptn += 1
            } else {
                gptu += 1
            }
        } else {
            if (traitPositive === 'true') {
                gntp += 1
            } else if (traitPositive === 'false') {
                gntn += 1
            } else {
                gntu += 1
            }
        }
    }
    return [gptp, gptn, gptu, gntp, gntn, gntu]
}


/**
 * Plot pie chart.
 *
 * @param {Element} targetElement Container element for plot.
 * @param {float|int} heightMultiplier The height will be scaled by this value
 * @param {[int]} values The values of the pie chart: ['g+t+', "g+t-", "g+t?", "g-t+", "g-t-", "g-t?"]
 * @param {string: string} colorDict Maps groups to colors, e.g. {"g+t+": "#1f78b4"}
 */
const plotPie = function ({targetElement, heightMultiplier, values, colorDict}) {
    console.log('plotPie:colorDict', colorDict)
    Plotly.purge(targetElement)

    const colNames = ['g+t+', "g+t-", "g+t?", "g-t+", "g-t-", "g-t?"]
    const [gptp, gptn, gptu, gntp, gntn, gntu] = values

    const data = [{
        values: [gptp, gptn, gptu, gntp, gntn, gntu],
        labels: colNames,
        type: 'pie',
        marker: {colors: colNames.map(col => colorDict[col])}
    }]
    const layout = calcLayout({}, heightMultiplier)
    Plotly.newPlot(targetElement, data, layout)
}


/**
 * Plot histogram.
 *
 * @param {Element} targetElement Container element for plot.
 * @param {float|int} heightMultiplier The height will be scaled by this value
 * @param {[[float]]} values The values of the pie chart: ['g+t+', "g+t-", "g+t?", "g-t+", "g-t-", "g-t?"]
 * @param {string: string} colorDict Maps groups to colors, e.g. {"g+t+": "#1f78b4"}
 */
const plotHistogram = function ({targetElement, heightMultiplier, values, colorDict}) {
    Plotly.purge(targetElement)

    const [gptp, gptn, gptu, gntp, gntn, gntu] = values

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
    }, heightMultiplier)

    Plotly.newPlot(targetElement, data, layout)
}


/**
 * geneClickFunctions: List of functions to be executed when a gene is clicked.
 * applyClick: runs the geneClickFunctions
 */
const geneClickFunctions = []
const applyClick = function (event) {
    const currentGene = event.target.textContent
    geneClickFunctions.forEach((fn) => fn(event, currentGene))
}


// Add click listeners to genes
Promise.all([tablePromise, coverageMatrixTablePromise])
    .then(([tableData, cmValues]) => {
        $('#coverage-matrix_wrapper thead th:not(:first):not(:last)').click((event) => {
            const gene = event.target.innerHTML
            if (!tableData.index.includes(gene)) return // ignore potential additional columns
            applyClick(event)
        })
        $('#result-tsv tbody tr td:first-child').click((event) => {
            applyClick(event)
        })
    })


// Plot first gene by default
Promise.all([tablePromise, coverageMatrixTablePromise, valuesPromise, configPromise])
    .then(([tableData, cmValues, valuesData, config]) => {
        const heightMultiplier = cmValues.isNumeric ? 0.5 : 1
        const targetElement = document.getElementById('pie')
        const targetNameElement = document.getElementById('pie-gene-name')

        const plotWrapper = function (event, currentGene) {
            targetNameElement.innerHTML = currentGene
            plotPie({
                targetElement: targetElement,
                heightMultiplier: heightMultiplier,
                values: calcPieValues(currentGene, cmValues.data, valuesData),
                colorDict: config['colors']
            })
        }

        // draw plot
        plotWrapper(null, tableData.index[0])

        // append to click event
        geneClickFunctions.push(plotWrapper)

        // add resize event
        let waitUntilResizeOverTimeout
        window.addEventListener('resize', () => {
            clearTimeout(waitUntilResizeOverTimeout)
            waitUntilResizeOverTimeout = setTimeout(plotWrapper, 300)
        })
    })


Promise.all([tablePromise, metaPromise, valuesPromise, coverageMatrixPromise, configPromise])
    .then(([tableData, metaData, valuesData, cmValues, config]) => {
        if (!valuesData.isNumeric) {
            console.log('no numeric data')
            return
        }

        const targetElement = document.getElementById('histogram')
        const targetNameElement = document.getElementById('histogram-gene-name')

        const plotWrapper = function (event, currentGene) {
            targetNameElement.innerHTML = currentGene
            plotHistogram({
                targetElement: targetElement,
                heightMultiplier: 0.5,
                values: calcHistValues(currentGene, cmValues.data, valuesData),
                colorDict: config['colors']
            })
        }

        // draw plot
        plotWrapper(null, tableData.index[0])

        // append to click event
        geneClickFunctions.push(plotWrapper)

        // add resize event
        let waitUntilResizeOverTimeout
        window.addEventListener('resize', () => {
            clearTimeout(waitUntilResizeOverTimeout)
            waitUntilResizeOverTimeout = setTimeout(plotWrapper, 300)
        })
    })


Promise.all([coverageMatrixPromise, valuesPromise, configPromise])
    .then(([cmValues, valuesData, config]) => {
        console.log('cmValues', cmValues)
        const concatString = config['concat-string'] ?? '+'
        console.log('valuesData', valuesData)
        // break if genes list
        if (!cmValues.isGeneList) return
        const genesLinks = config['genes-links']
        // break if no links
        if (genesLinks.length === 0) return

        /**
         * Generate content of popover.
         */
        const genesPopoverContent = function () {
            const orthoGene = this.textContent
            const orthoGeneCoord = cmValues.rawGenes.indexOf(orthoGene)
            const orthoGeneSanitized = cmValues.sanitizedGenes[orthoGeneCoord]
            const context = {
                'orthogene': orthoGene,
                'have-isolates': [],
                'lack-isolates': [],
                'all-genes': [],
                'positive-genes': [],
                'negative-genes': [],
                'unclear-genes': []
            }
            for (const isolateData of cmValues.geneList) {
                const isolate = isolateData['Isolate']
                const genes = isolateData[orthoGeneSanitized]
                const traitPositive = valuesData.isolateToTrait[isolate]

                if (genes.length === 0) {
                    context['lack-isolates'].push(isolate)
                    continue
                }
                context['all-genes'].push(...genes)
                context['have-isolates'].push(isolate)

                if (traitPositive === 'true') {
                    context['positive-genes'].push(...genes)
                } else if (traitPositive === 'false') {
                    context['negative-genes'].push(...genes)
                } else {
                    context['unclear-genes'].push(...genes)
                }
            }

            // concat arrays with config.json['concat-string']
            for (const [key, value] of Object.entries(context)) {
                if (Array.isArray(value)) {
                    context[key] = value.join(concatString)
                }
            }

            let html = '<div class="list-group">'

            for (const [linkName, linkTemplate] of Object.entries(genesLinks)) {
                const link = linkTemplate.format(context)
                html += `<a href="${link}" class="list-group-item list-group-item-action">${linkName}<\a>`
            }
            html += '</div>'

            return html
        }

        const genesPopover = function (event, currentGene) {
            console.log('0', event)
            console.log('1', currentGene)
            console.log('2', cmValues)

            event.target.classList.add('has-popover')

            new bootstrap.Popover(event.target, {
                container: 'body',
                trigger: 'click',
                html: true,
                title: currentGene,
                content: genesPopoverContent
            }).show()

            popovers.push(event.target)
        }

        // append to click event
        geneClickFunctions.push(genesPopover)


    })


// load newick file (tree)
const newickPromise = fetch('tree.nwk')
    .then(response => response.text())


const getSize = function (element) {
    const computedStyle = getComputedStyle(element)
    return element.clientWidth - parseFloat(computedStyle.paddingLeft) - parseFloat(computedStyle.paddingRight)
}
// draw tree
const treePromise = Promise.all([documentReadyPromise, newickPromise, valuesPromise, configPromise])
    .then(([_, newickTree, valuesData, config]) => {
        // https://www.phylocanvas.gl/examples/metadata-blocks.html
        // https://www.phylocanvas.gl/docs/methods.html
        const parentDiv = document.getElementById('tree-and-tables')
        const treeContainer = document.querySelector("#tree")
        let currentWidth = getSize(parentDiv)
        const tree = new phylocanvas.PhylocanvasGL(
            treeContainer, {
                type: phylocanvas.TreeTypes.Hierarchical,
                blocks: valuesData.isNumeric ? ['gene', 'trait', 'value'] : ['gene', 'trait'],
                blockLength: 24, // block size in pixels
                blockPadding: 8, // the gap size between blocks in pixels
                alignLabels: true,
                showLabels: true,
                showLeafLabels: true,
                size: {width: currentWidth, height: config['tree-height']},
                source: newickTree,
            }
        )

        const redrawWrapper = function () {
            const newWidth = getSize(parentDiv)
            if (newWidth !== currentWidth) {
                tree.resize(newWidth, config['tree-height'])
                currentWidth = newWidth
            }
        }
        let waitUntilResizeOverTimeout
        window.addEventListener('resize', () => {
            clearTimeout(waitUntilResizeOverTimeout)
            waitUntilResizeOverTimeout = setTimeout(redrawWrapper, 300)
        })

        return tree
    })


/**
 * Update the colorful rows of metadata below the phylogenetic tree
 *
 * @param {Object} tree phylocanvas.gl object
 * @param {string} gene orthogene name (unsanitized)
 * @param {[]} coverageMatrix Gene count table
 * @param {[]} valuesData Boolean trait table
 * @param {string: string} colorDict Maps groups to colors, e.g. {"g+t+": "#1f78b4"}
 */
const updateTreeMetadata = function (tree, gene, coverageMatrix, valuesData, colorDict) {
    const sanitizedGene = sanitizeKey(gene)

    // calculate metadata
    const metadata = {}
    for (const isolateData of coverageMatrix) {
        const isolate = isolateData['Isolate']
        const genePositive = isolateData[sanitizedGene].toString() !== '0'
        const traitPositive = valuesData.isolateToTrait[isolate]
        // const traitValue = valuesData.isolateToValue[isolate]

        metadata[isolate] = {trait: {}, gene: {}, value: {}}

        if (genePositive) {
            metadata[isolate]['gene']['colour'] = colorDict['g+']
        } else {
            metadata[isolate]['gene']['colour'] = colorDict['g-']
        }

        if (traitPositive === 'true') {
            metadata[isolate]['trait']['colour'] = colorDict['t+']
        } else if (traitPositive === 'false') {
            metadata[isolate]['trait']['colour'] = colorDict['t-']
        } else {
            metadata[isolate]['trait']['colour'] = colorDict['t?']
        }

        if (valuesData.isNumeric) {
            metadata[isolate]['value']['colour'] = valuesData.colorScale(valuesData.isolateToNormValue[isolate])
        }
    }

    // update tree
    tree.setProps({
        metadata: metadata
    })

    return metadata
}


Promise.all([treePromise, valuesPromise, coverageMatrixPromise, tablePromise, configPromise])
    .then(([tree, valuesData, cmValues, tableData, config]) => {
        const targetNameElement = document.getElementById('tree-gene-name')

        const plotWrapper = function (event, currentGene) {
            targetNameElement.innerHTML = currentGene
            updateTreeMetadata(tree, currentGene, cmValues.data, valuesData, config['colors'])
        }

        // draw plot
        plotWrapper(null, tableData.index[0])

        // append to click event
        geneClickFunctions.push(plotWrapper)
    })
