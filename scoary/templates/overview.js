"use strict"
for (const attr of ['table', 'thead', 'tbody', 'tr', 'td', 'th']) {
    bootstrap.Tooltip.Default.allowList[attr] = []
}

const svgContainer = document.getElementById('svg-container')
let svgDocument
let overviewDf
let patches


const loadSvgPromise = fetch('overview_plot.svg', {method: 'GET', headers: {}}).then(function (response) {
    return response.text()
}).then((data) => {
    svgContainer.innerHTML = data
    svgDocument = svgContainer.firstElementChild
    overviewDf = JSON.parse(svgContainer.childNodes[4].data)
    patches = svgDocument.getElementById('clickable-patches')
}).catch((err) => {
    alert('Failed to load overview_plot.svg')
    console.log(this, err)
})


function createTooltipTitle() {
    const elementIndex = getElementIndex(this)
    return overviewDf['index'][elementIndex]
}

function createTooltipContent() {
    const elementIndex = getElementIndex(this)
    const traitName = overviewDf['index'][elementIndex]

    let html = '<table id="popoverTable" class="table"><tbody>'

    overviewDf["columns"].forEach((x, i) => {
        const title = x.replace(/^(min_)/, '')
        const value = overviewDf['data'][elementIndex][i]
        console.log(title, value)
        html += `
                <tr>
                    <th scope="row">${title}</th>
                    <td>${value.toPrecision(2)}</td>
                </tr>`
    })

    html += '</tbody></table>'
    return html
}

// activate clickable boxes popovers
loadSvgPromise.then(function (result) {
        Array.from(patches.children).forEach((patch, index) => {
            new bootstrap.Popover(patch, {
                container: 'body',
                trigger: 'hover',
                html: true,
                title: createTooltipTitle,
                content: createTooltipContent
            })
        })
    }
)

// // activate clickable boxes
loadSvgPromise.then(function (result) {
        Array.from(patches.children).forEach((patch, index) => {
            patch.addEventListener("click", redirect)
        })
    }
)

// redirect to traits page
function redirect(event) {
    const elementIndex = getElementIndex(this)
    const traitName = overviewDf['index'][elementIndex]
    console.log(`FORWARDING TO: ${traitName}`)
    window.location.href = `./trait.html?trait=${traitName}`;
}

function getElementIndex(element) {
    return [...element.parentNode.children].indexOf(element);
}


// function Counter(array) {
//     let count = {}
//     array.forEach(val => count[val] = (count[val] || 0) + 1)
//     return count
// }
//
// function getMsCpd(el) {
//     return el.firstElementChild.innerHTML.split(' -->', 1)[0].split('<!-- ')[1].split(':')
// }
//
// function highlightCpd(cpdName) {
//     console.log(cpdName)
//     for (var el of document.getElementById('matplotlib.axis_2').children) {
//         const currCpdName = getMsCpd(el)[1]
//         if (currCpdName === cpdName) {
//             el.style['stroke'] = 'yellow'
//             el.style['stroke-width'] = '1px'
//         } else {
//             el.style['stroke'] = ''
//             el.style['stroke-width'] = ''
//         }
//     }
// }
//
// function clickOnCpd(e) {
//     highlightCpd(getMsCpd(e.target.parentElement.parentElement)[1])
// }
//
// const cpds = Array.prototype.slice.call(document.getElementById('matplotlib.axis_2').children).map(el => getMsCpd(el)[1])
// const cpdCounts = Counter(cpds)
//
// for (var el of document.getElementById('matplotlib.axis_2').children) {
//     const data = getMsCpd(el)
//     const cpdName = data[1]
//     const dataSet = data[0]
//
//     if (cpdCounts[cpdName] > 1) {
//         el.style['fill'] = 'darkred'
//     }
//
//     el.addEventListener('contextmenu', function (e) {
//         clickOnCpd(e)
//     })
//
//     bbox = el.getBBox()
//     let a = document.createElementNS("http://www.w3.org/2000/svg", "a")
//     a.setAttribute('href', `${dataSet}/gm_prob/${cpdName}-genes.html`)
//
//     let newRect = document.createElementNS("http://www.w3.org/2000/svg", "rect")
//     newRect.setAttribute("x", bbox.x)
//     newRect.setAttribute("y", bbox.y)
//     newRect.setAttribute("width", bbox.width)
//     newRect.setAttribute("height", bbox.height)
//     newRect.setAttribute("fill", "transparent")
//     newRect.setAttribute("data-compound", cpdName)
//     newRect.addEventListener('click', function () {
//         console.log(cpdName)
//     })
//
//     a.append(newRect)
//     el.append(a)  // in front
// }
