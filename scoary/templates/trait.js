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
document.title = `${trait} (Scoary2)`


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
documentReadyPromise.then(() => {
    document.querySelector('#trait-name').textContent = trait
    // set the download links
    document.getElementById('download-genes').href = `traits/${trait}/result.tsv`
    document.getElementById('download-coverage-matrix').href = `traits/${trait}/coverage-matrix.tsv`
    document.getElementById('download-values').href = `traits/${trait}/values.tsv`
})


/**
 * Promise for config.json. Returns the dictionary object.
 */
const configPromise = fetch(`app/config.json`)
    .then(response => response.json())
    .then(config => {
        if (config['table-config']['sanitize-genes']) {
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
 * Promise for tree.nwk. Returns the dictionary object.
 */
const newickPromise = fetch('tree.nwk')
    .then(response => response.text())


/**
 * Show metadata.
 */
metaPromise.then(metaData => {
    // save default metadata
    document.getElementById('metadata-content').textContent = JSON.stringify(metaData, null, '\t')


    if (metaData.hasOwnProperty('info')) {
        // add trait info

        let html = '<table id="popoverTable" class="table"><tbody>'

        for (const [key, value] of Object.entries(metaData.info)) {
            html += `
                <tr>
                    <th scope="row"><code>${key}</code></th>
                    <td><code>${value}</code></td>
                </tr>`
        }

        document.getElementById('trait-content').innerHTML = html
        // show element
        document.getElementById('trait-metadata').hidden = false
        // do not show metaData.info twice
        delete metaData.info

        html += '</tbody></table>'

    }

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
        if (!valuesData.isNumeric) {
            document.getElementById('histogram-container').remove()
            document.getElementById('download-values').style.display = 'none'
        }
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
            valuesData.colorScale = chroma.scale(config['tree-config']['color-scale'])
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
        const tableConfig = config['table-config']
        const floatCols = tableConfig['float-cols']
        const renderFloat = (data, type, row) => parseFloat(data).toPrecision(3)
        tableData.columns = tableData.meta.fields.map(col => {
            const column = {'data': col, 'label': col, 'title': col}
            if (floatCols.includes(col)) column['render'] = renderFloat
            return column
        })
        tableData.index = tableData['data'].map(col => col['Gene'])
        tableData.orderCol = tableData.meta.fields.indexOf(
            tableData.meta.fields.includes('fq*ep') ? 'fq*ep' : 'fisher_q')
        tableData.order = 'asc'
        tableData.hiddenCols = tableConfig['default-hidden-cols'].map(col => tableData.meta.fields.indexOf(col))

        console.log('tableData', tableData)

        return tableData
    })
})


/**
 * Load Genes table.
 */
const genesTablePromise = Promise.all([tablePromise, documentReadyPromise]).then(([tableData, _]) => {
    console.log(tableData.data)
    console.log(tableData.columns)
    console.log(tableData.hiddenCols)
    tableData.table = $('#result-tsv').DataTable({
        data: tableData.data,
        columns: tableData.columns,
        paging: false,
        responsive: true,
        order: [[tableData.orderCol, tableData.order]],
        columnDefs: [{
            "targets": tableData.hiddenCols, "visible": false, "searchable": false
        }]
    })
    return tableData
}).then((tableData) => {
    // load buttons to toggle columns
    const target = document.getElementById('vis-toggler')
    let htmls = []
    tableData.meta.fields.forEach((col, i) => {
        const active = tableData.hiddenCols.includes(i) ? '' : 'active'
        htmls.push(`<button class="toggle-vis btn btn-sm btn-primary ${active}" data-column="${i}">${col}</button>`)
    })
    target.innerHTML = htmls.join(' ')
    $('.toggle-vis').on('click', function (event) {
        const column = tableData.table.column($(this).attr('data-column'))
        if (column.visible()) {
            column.visible(false)
            event.target.classList.remove('active')
        } else {
            column.visible(true)
            event.target.classList.add('active')
        }
    })
    return tableData
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
    console.info('no isolate info')
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
const hidePopovers = () => {
    while (popovers.length > 0) {
        const popover = popovers.pop()
        bootstrap.Popover.getInstance(popover).hide()
    }
}
$(document).on("click", (event) => {
    // ignore clicks on .popover
    if (!Boolean(event.target.closest('.popover, table'))) {
        hidePopovers()
    }
})


/**
 * Load the coverage matrix table, make it interactive.
 */

const coverageMatrixTablePromise = Promise.all(
    [coverageMatrixPromise, configPromise]
).then(([cmValues, config]) => {
    // Load DataTable
    const generateColumnDef = (col, i) => {
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

    if (!cmValues.isGeneList) return cmValues

    /**
     * Make genes clickable depending on config.json.
     */
    const createGeneElement = (gene) => {
        let geneLink = config['link-config']['single-gene']
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
    const allGeneCells = document.querySelectorAll('#coverage-matrix tbody')

    document.getElementById('coverage-matrix').addEventListener('click', (event) => {
        if (!event.target.classList.contains('gene-count')) return
        let instance = bootstrap.Popover.getInstance(event.target)
        if (instance === null) {
            instance = new bootstrap.Popover(event.target, {
                container: '#container-left',
                trigger: 'manual',
                html: true,
                title: genePopoverTitle,
                content: genePopoverContent
            })
        }
        instance.show()
    })


    return cmValues
})


/**
 * Plot pie chart.
 *
 * @param {Element} targetElement Container element for plot.
 * @param {[int]} values The values of the pie chart: ['g+t+', "g+t-", "g+t?", "g-t+", "g-t-", "g-t?"]
 * @param {string: string} colorDict Maps groups to colors, e.g. {"g+t+": "#1f78b4"}
 */
const plotPie = ({targetElement, values, colorDict}) => {
    Plotly.purge(targetElement)

    const colNames = ['g+t+', "g+t-", "g+t?", "g-t+", "g-t-", "g-t?"]
    const [gptp, gptn, gptu, gntp, gntn, gntu] = values

    const data = [{
        values: [gptp, gptn, gptu, gntp, gntn, gntu],
        labels: colNames,
        type: 'pie',
        marker: {colors: colNames.map(col => colorDict[col])}
    }]
    Plotly.newPlot(targetElement, data)
}


/**
 * Plot histogram.
 *
 * @param {Element} targetElement Container element for plot.
 * @param {[[float]]} values The values of the pie chart: ['g+t+', "g+t-", "g+t?", "g-t+", "g-t-", "g-t?"]
 * @param {string: string} colorDict Maps groups to colors, e.g. {"g+t+": "#1f78b4"}
 */
const plotHistogram = ({targetElement, values, colorDict}) => {
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

    const layout = {barmode: "stack", xaxis: {title: "numerical value"}, yaxis: {title: "count"}}

    Plotly.newPlot(targetElement, data, layout)
}


/**
 * geneClickFunctions: List of functions to be executed when a gene is clicked.
 * applyClick: runs the geneClickFunctions
 */
const geneClickFunctions = []
const applyClick = (event) => {
    const currentGene = event.target.textContent
    Promise.all([coverageMatrixDataPromise, valuesPromise])
        .then(([cmValues, valuesData]) => {
            const orthogeneData = calcOrthogeneData(currentGene, cmValues.data, valuesData)
            geneClickFunctions.forEach((fn) => fn(event, currentGene, orthogeneData))
        })
}


// Add click listeners to genes
Promise.all([tablePromise, coverageMatrixTablePromise])
    .then(([tableData, _]) => {
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
        const targetElement = document.getElementById('pie')
        const targetNameElement = document.getElementById('pie-gene-name')

        const plotWrapper = (event, currentGene, orthogeneData) => {
            targetNameElement.innerHTML = currentGene
            plotPie({
                targetElement: targetElement,
                values: orthogeneData.map(x => x.length),
                colorDict: config['colors']
            })
        }

        // draw plot
        const currentGene = tableData.index[0]
        plotWrapper(null, currentGene, calcOrthogeneData(currentGene, cmValues.data, valuesData))

        // append to click event
        geneClickFunctions.push(plotWrapper)

        // add resize event
        let waitUntilResizeOverTimeout
        window.addEventListener('resize', () => {
            clearTimeout(waitUntilResizeOverTimeout)
            waitUntilResizeOverTimeout = setTimeout(() => Plotly.Plots.resize(targetElement), 300)
        })
    })


Promise.all([tablePromise, metaPromise, valuesPromise, coverageMatrixPromise, configPromise])
    .then(([tableData, metaData, valuesData, cmValues, config]) => {
        if (!valuesData.isNumeric) {
            console.info('no numeric data')
            return
        }

        const targetElement = document.getElementById('histogram')
        const targetNameElement = document.getElementById('histogram-gene-name')

        const plotWrapper = (event, currentGene, orthogeneData) => {
            targetNameElement.innerHTML = currentGene
            plotHistogram({
                targetElement: targetElement,
                values: orthogeneData,
                colorDict: config['colors']
            })
        }

        // draw plot
        const currentGene = tableData.index[0]
        plotWrapper(null, currentGene, calcOrthogeneData(currentGene, cmValues.data, valuesData))

        // append to click event
        geneClickFunctions.push(plotWrapper)

        // add resize event
        let waitUntilResizeOverTimeout
        window.addEventListener('resize', () => {
            clearTimeout(waitUntilResizeOverTimeout)
            waitUntilResizeOverTimeout = setTimeout(() => Plotly.Plots.resize(targetElement), 300)
        })
    })


Promise.all([coverageMatrixPromise, valuesPromise, configPromise])
    .then(([cmValues, valuesData, config]) => {
        const concatString = config['link-config']['concat-string'] ?? '+'
        // break if genes list
        if (!cmValues.isGeneList) return
        const genesLinks = config['link-config']['many-genes']
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

        const genesPopover = (event, currentGene) => {
            let instance = bootstrap.Popover.getInstance(event.target)
            if (instance === null) {
                instance = new bootstrap.Popover(event.target, {
                    container: '#container-left',
                    trigger: 'manual',
                    html: true,
                    title: currentGene,
                    content: genesPopoverContent
                })
            }
            instance.show()

            popovers.push(event.target)
        }

        // append to click event
        geneClickFunctions.push(genesPopover)
    })


const getInnerWidth = (element) => {
    const computedStyle = getComputedStyle(element)
    return element.clientWidth - parseFloat(computedStyle.paddingLeft) - parseFloat(computedStyle.paddingRight)
}


// draw tree
const treePromise = Promise.all([documentReadyPromise, newickPromise, valuesPromise, configPromise])
    .then(([_, newickTree, valuesData, config]) => {
        // https://www.phylocanvas.gl/examples/metadata-blocks.html
        // https://www.phylocanvas.gl/docs/methods.html
        const treeContainer = document.querySelector("#tree")
        let currentWidth = getInnerWidth(treeContainer.parentElement)
        const tree = new phylocanvas.PhylocanvasGL(
            treeContainer, {
                type: phylocanvas.TreeTypes[config['tree-config']['type']],
                blocks: valuesData.isNumeric ? ['gene', 'trait', 'value'] : ['gene', 'trait'],
                blockLength: 24, // block size in pixels
                blockPadding: 8, // the gap size between blocks in pixels
                alignLabels: true,
                showLabels: true,
                showLeafLabels: true,
                interactive: {zoom: false, pan: false, tooltip: true, highlight: true},
                zoom: false,
                size: {width: currentWidth, height: config['tree-config']['height']},
                source: newickTree,
            }
        )

        const redrawWrapper = () => {
            const newWidth = getInnerWidth(treeContainer.parentElement)
            if (newWidth !== currentWidth) {
                tree.resize(newWidth, config['tree-config']['height'])
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
 * @param {[]} coverageMatrix Gene count table
 * @param {[]} valuesData Boolean trait table
 * @param {string: string} colorDict Maps groups to colors, e.g. {"g+t+": "#1f78b4"}
 * @returns {string: Array} metadata Dictionary of metadata
 */
const loadTreeMetadata = (coverageMatrix, valuesData, config) => {
    const metaConfig = config['tree-config']['metadata-bars']

    // calculate metadata
    const metadata = {}
    for (const isolateData of coverageMatrix) {
        const isolate = isolateData['Isolate']
        const traitPositive = valuesData.isolateToTrait[isolate]

        let isolateConfig = {}
        if (traitPositive === 'true') {
            isolateConfig['trait'] = metaConfig['t+']
        } else if (traitPositive === 'false') {
            isolateConfig['trait'] = metaConfig['t-']
        } else {
            isolateConfig['trait'] = metaConfig['t?']
        }

        if (valuesData.isNumeric) {
            isolateConfig['value'] = {
                colour: valuesData.colorScale(valuesData.isolateToNormValue[isolate]),
                label: parseFloat(valuesData.isolateToValue[isolate]).toPrecision(3)
            }
        }

        metadata[isolate] = isolateConfig
    }

    return metadata
}

// https://www.phylocanvas.gl/examples/styling-nodes.html
const updateTreeGeneColors = (tree, gene, baseMetadata, coverageMatrix, valuesData, config) => {
    const sanitizedGene = sanitizeKey(gene)
    const metaConfig = config['tree-config']['metadata-bars']
    const leafConfig = config['tree-config']['leaf-nodes']

    // calculate styles
    const styles = {}
    const metadata = baseMetadata
    for (const isolateData of coverageMatrix) {
        const isolate = isolateData['Isolate']
        const genePositive = isolateData[sanitizedGene].toString() !== '0'

        styles[isolate] = {}

        if (genePositive) {
            styles[isolate] = leafConfig['g+']
            metadata[isolate]['gene'] = metaConfig['g+']
        } else {
            styles[isolate] = leafConfig['g-']
            metadata[isolate]['gene'] = metaConfig['g-']
        }

    }

    // update tree
    tree.setProps({
        styles: styles,
        metadata: metadata
    })

    return styles
}

/**
 * For a gene, calculate the numeric values per category. (Raw data for histogram.)
 *
 * @param {string} gene Name of the orthgene.
 * @param {[]} coverageMatrix Gene count table
 * @param {[]} valuesData Boolean trait table
 */
const calcOrthogeneData = (gene, coverageMatrix, valuesData) => {
    const sanitizedGene = sanitizeKey(gene)
    const [gptp, gptn, gptu, gntp, gntn, gntu] = [[], [], [], [], [], []]
    for (const isolateData of coverageMatrix) {
        const isolate = isolateData['Isolate']
        const genePositive = isolateData[sanitizedGene].toString() !== '0'
        const traitPositive = valuesData.isolateToTrait[isolate]
        const traitValue = valuesData.isNumeric ? valuesData.isolateToValue[isolate] : null

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


Promise.all([treePromise, valuesPromise, coverageMatrixPromise, tablePromise, configPromise])
    .then(([tree, valuesData, cmValues, tableData, config]) => {
        const targetNameElement = document.getElementById('tree-gene-name')

        // load colorful bars below tree
        const baseMetadata = loadTreeMetadata(cmValues.data, valuesData, config)

        const plotWrapper = (event, currentGene) => {
            targetNameElement.innerHTML = currentGene
            updateTreeGeneColors(tree, currentGene, baseMetadata, cmValues.data, valuesData, config)
        }

        // draw plot
        plotWrapper(null, tableData.index[0])

        // append to click event
        geneClickFunctions.push(plotWrapper)
    })
