"use strict"


for (const attr of ['table', 'thead', 'tbody', 'tr', 'td', 'th']) {
    bootstrap.Tooltip.Default.allowList[attr] = []
}

const svgContainer = document.getElementById('svg-container')
const traitInfoDiv = document.getElementById('trait-info')
let svgDocument, overviewDf, patches


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


const isFloat = (n) => {
    n = parseFloat(n)
    if (Number.isNaN(n)) return false
    return Number(n) === n && n % 1 !== 0
}

function createTooltipContent(elementIndex) {
    const traitData = overviewDf.data[elementIndex]
    const traitName = overviewDf['index'][elementIndex]

    let html = `<h2><a href="trait.html?trait=${traitName}">${traitName} 
<button type="button" class="btn btn-info"><svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-box-arrow-up-right" viewBox="0 0 16 16">
  <path fill-rule="evenodd" d="M8.636 3.5a.5.5 0 0 0-.5-.5H1.5A1.5 1.5 0 0 0 0 4.5v10A1.5 1.5 0 0 0 1.5 16h10a1.5 1.5 0 0 0 1.5-1.5V7.864a.5.5 0 0 0-1 0V14.5a.5.5 0 0 1-.5.5h-10a.5.5 0 0 1-.5-.5v-10a.5.5 0 0 1 .5-.5h6.636a.5.5 0 0 0 .5-.5z"></path>
  <path fill-rule="evenodd" d="M16 .5a.5.5 0 0 0-.5-.5h-5a.5.5 0 0 0 0 1h3.793L6.146 9.146a.5.5 0 1 0 .708.708L15 1.707V5.5a.5.5 0 0 0 1 0v-5z"></path>
</svg></button>
</a></h2>`

    html += '<table id="popoverTable" class="table"><tbody>'

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
    const patchesElement = document.getElementById('clickable-patches')

    patchesElement.addEventListener('mouseover', function (event) {
        traitInfoDiv.innerHTML = createTooltipContent(event.target.getAttribute('index'))
    })
    patchesElement.addEventListener('auxclick', openTraitsTab)
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

let ax4, slimSelect
const selected = []

// toggle trait yellow
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

const goToTrait = () => {
    const selectId = parseInt(slimSelect.selected())
    const target = ax4.children[selectId]
    target.scrollIntoView()
}

// // activate clickable boxes
Promise.all([loadSvgPromise, overviewPromise]
).then(([result, __]) => {
    const copyButton = document.getElementById('copy-button')
    copyButton.addEventListener('click', copySelectedToClipboard)

    const deselectButton = document.getElementById('deselect-button')
    deselectButton.addEventListener('click', deselectAll)


    ax4 = document.getElementById('clickable-patches')

    slimSelect = new SlimSelect({
        select: '#slim-select',
        placeholder: 'Click here to select traits',
        data: overviewDf.index.map((trait, index) => {
            return {text: trait, value: index.toString()}
        }),
        onChange: function () {
            const selectId = parseInt(slimSelect.selected())
            toggleTrait(selectId)
            console.log(selectId)
            traitInfoDiv.innerHTML = createTooltipContent(selectId)
        },
    })

    const goToButton = document.getElementById('goto-button')
    goToButton.addEventListener('click', goToTrait)
})
