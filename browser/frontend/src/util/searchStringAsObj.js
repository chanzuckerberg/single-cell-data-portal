const searchStringAsObj = search => {
  let obj = null
  if (search) {
    obj = JSON.parse(
      '{"' +
        decodeURI(search)
          .replace(/"/g, '\\"')
          .replace(/&/g, '","')
          .replace(/=/g, '":"') +
        '"}'
    )
  }
  return obj
}

export default searchStringAsObj
