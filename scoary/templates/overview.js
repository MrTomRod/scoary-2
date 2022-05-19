"use strict"


for (const attr of ['table', 'thead', 'tbody', 'tr', 'td', 'th']) {
    bootstrap.Tooltip.Default.allowList[attr] = []
}

const svgContainer = document.getElementById('svg-container')
let svgDocument
let overviewDf
let patches


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
const overviewPromise = Papa.execPromise('summary.tsv', {
    header: true, download: true, skipEmptyLines: true, delimiter: '\t', newline: '\n',
}).then(overviewData => {
    // Preprocess data
    overviewData.index = overviewData['data'].map(col => col['Trait'])
    overviewDf = overviewData
    return overviewData
}).catch(() => {
    console.info('no isolate info')
    return 'no isolate info'
})

overviewPromise.then((overviewData) => {
    if (overviewData.data.length == 1) {
        const traitName = overviewData.data[0]['Trait']
        const url = `./trait.html?trait=${traitName}`
        console.info(`FORWARDING TO: ${url}`)
        window.location.href = url
    }
})

const loadSvgPromise = fetch('overview_plot.svg', {method: 'GET', headers: {}}).then(function (response) {
    return response.text()
}).then((data) => {
    svgContainer.innerHTML = data
    svgDocument = svgContainer.firstElementChild
    patches = svgDocument.getElementById('clickable-patches')
}).catch((err) => {
    alert('Failed to load overview_plot.svg')
    console.warn(this, err)
})

function createTooltipTitle() {
    const elementIndex = this.getAttribute('index')
    const traitName = overviewDf['index'][elementIndex]
    return `<a href="trait.html?trait=${traitName}">${traitName}</a>`
}

const isFloat = (n) => {
    n = parseFloat(n)
    if (Number.isNaN(n)) return false
    return Number(n) === n && n % 1 !== 0
}

function createTooltipContent() {
    const elementIndex = this.getAttribute('index')
    const traitData = overviewDf.data[elementIndex]

    let html = '<table id="popoverTable" class="table"><tbody>'

    overviewDf.meta.fields.forEach((key, i) => {
        if (key === 'Trait') return
        let value = traitData[key]
        if (value === '') return
        key = key.replace(/^(min_)/, '')
        if (isFloat(value)) {
            value = parseFloat(value).toPrecision(2)
        }

        html += `
            <tr>
                <th scope="row">${key}</th>
                <td>${value}</td>
            </tr>`
    })

    html += '</tbody></table>'
    return html
}

Promise.all([loadSvgPromise, overviewPromise]).then(() => {
    const openPopovers = []
    const patchesElement = document.getElementById('clickable-patches')

    const hideAllPopovers = () => {
        while (openPopovers.length) {
            const target = openPopovers.pop()
            bootstrap.Popover.getInstance(target).hide()
        }
    }
    patchesElement.addEventListener('mouseover', function (event) {
        let instance = bootstrap.Popover.getInstance(event.target)
        if (instance === null) {
            instance = new bootstrap.Popover(event.target, {
                container: 'body',
                trigger: 'manual',
                placement: 'right',
                html: true,
                title: createTooltipTitle,
                content: createTooltipContent,
            })
        }

        instance.show()
        hideAllPopovers()
        openPopovers.push(event.target)
    })
    patchesElement.addEventListener('auxclick', openTraitsTab)
    document.addEventListener('click', function (event) {
        const isPopoverOrSvg = event.target.closest('.popover, svg') !== null
        if (isPopoverOrSvg) return
        hideAllPopovers()
    })
    Array.from(patchesElement.children).forEach((patch, index) => {
        patch.setAttribute('index', index)
    })
    patchesElement.addEventListener('click', toggleTrait)
})


// Open traits page in new tab
function openTraitsTab(event) {
    const traitId = event.target.getAttribute('index')
    const traitName = overviewDf.index[traitId]
    const url = `./trait.html?trait=${traitName}`
    console.info(`FORWARDING TO: ${url}`)
    // window.location.href = url
    window.open(url, '_blank').focus()
}

let ax4
const selected = []

const toggleTrait = (eventOrId) => {
    let target, index
    if (Number.isInteger(eventOrId)) {
        index = eventOrId
        target = ax4.children[eventOrId]
    } else {
        target = eventOrId.target
        index = target.getAttribute('index')
    }
    if (selected.includes(index)) {
        for (var i = 0; i < selected.length; i++) {
            if (selected[i] === index) selected.splice(i, 1)
        }
        target.style['fill'] = 'rgb(255, 255, 255)'
        target.style['fill-opacity'] = '0.000001'
    } else {
        selected.push(index)
        target.style['fill'] = 'yellow'
        target.style['fill-opacity'] = '0.2'
    }
}

const copySelectedToClipboard = () => {
    const content = selected.map(index => overviewDf.index[index]).join('\t')
    console.info('Copying to clipboard:', content)
    copyToClipboard(content)
}

function copyToClipboard(text) {
    // navigator clipboard api needs a secure context (https)
    if (navigator.clipboard && window.isSecureContext) {
        // navigator clipboard api method'
        return navigator.clipboard.writeText(text);
    } else {
        // text area method
        let textArea = document.createElement("textarea")
        textArea.value = text
        textArea.style.position = "fixed"
        textArea.style.left = "-999999px"
        textArea.style.top = "-999999px"
        document.body.appendChild(textArea)
        textArea.focus()
        textArea.select()
        return new Promise((res, rej) => {
            document.execCommand('copy') ? res() : rej()
            textArea.remove()
        })
    }
}

const deselectAll = () => {
    while (selected.length) {
        const target = ax4.children[selected.pop()]
        target.style['fill'] = 'rgb(255, 255, 255)'
        target.style['fill-opacity'] = '0.000001'
    }
}

// // activate clickable boxes
Promise.all([loadSvgPromise, overviewPromise]
).then(([result, __]) => {
    const copyButton = document.getElementById('copy-button')
    copyButton.addEventListener('click', copySelectedToClipboard)

    const deselectButton = document.getElementById('deselect-button')
    deselectButton.addEventListener('click', deselectAll)


    ax4 = document.getElementById('clickable-patches')

    const slimSelect = new SlimSelect({
        select: '#slim-select',
        placeholder: 'Click here to select traits',
        data: overviewDf.index.map((trait, index) => {
            return {text: trait, value: index.toString()}
        }),
        onChange: function () {
            const selectId = parseInt(slimSelect.selected())
            toggleTrait(selectId)
        },
    })
})
