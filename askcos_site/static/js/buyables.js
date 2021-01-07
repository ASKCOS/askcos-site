var app = new Vue({
    el: '#app',
    data: {
        buyables: [],
        uploadFile: '',
        searchSmilesQuery: '',
        searchSourceQuery: '',
        searchRegex: false,
        searchLimit: 100,
        canonSmiles: true,
        allowOverwrite: true,
        showAddModal: false,
        showUploadModal: false,
        addBuyableSmiles: '',
        addBuyablePrice: 1,
        addBuyableSource: '',
        uploadFileFormat: 'json'
    },
    mounted: function() {
        var urlParams = new URLSearchParams(window.location.search)
        let query = urlParams.get('q')
        if (!!query) {
            this.searchSmilesQuery = query
            this.search()
        }
    },
    methods: {
        search: function() {
            showLoader()
            fetch('/api/buyables/search?q='+encodeURIComponent(this.searchSmilesQuery)+'&source='+this.searchSourceQuery+'&regex='+this.searchRegex+'&limit='+this.searchLimit+'&canonicalize='+this.canonSmiles)
            .then(resp => resp.json())
            .then(json => {
                console.log(json['buyables']);
                this.buyables = json['buyables'];
            })
            .finally(() => hideLoader())
        },
        handleFileUpload: function() {
            this.uploadFile = this.$refs.file.files[0]
            if (this.uploadFile.name.endsWith('.json')) {
                this.uploadFileFormat = 'json'
            }
            if (this.uploadFile.name.endsWith('.csv')) {
                this.uploadFileFormat = 'csv'
            }
        },
        handleUploadSubmit: function() {
            showLoader()
            if (this.uploadFile == '') {
                alert('Please select a file to upload')
                hideLoader()
                return
            }
            let formData = new FormData()
            formData.append('file', this.uploadFile)
            formData.append('format', this.uploadFileFormat)
            formData.append('allowOverwrite', this.allowOverwrite)
            fetch('/api/buyables/upload', 
                {
                    'method': 'POST',
                    headers: {
                        'X-CSRFToken': getCookie('csrftoken')
                    },
                    body: formData
                }
            )
            .then(resp => resp.json())
            .then(json => {
                if (json.error) {
                    alert(json.error)
                    hideLoader()
                    return
                }
                alert('Out of '+json.total+' entries, successfully added '+json.added_count+', updated '+json.updated_count+', and skipped '+json.duplicate_count+' duplicates. Only adding (up to) '+2*this.searchLimit+' documents to the list below')
                if (json.added.length > 0) {
                    this.buyables.unshift(...json.added)
                }
                if (json.updated.length > 0) {
                    for (updated of json.updated) {
                        let inList = false
                        for (buyable of this.buyables) {
                            if (buyable._id == updated._id) {
                                inList = true
                                buyable.ppg = updated.ppg
                                buyable.source = updated.source
                                break
                            }
                        }
                        if (!inList) {
                            this.buyables.unshift(updated)
                        }
                    }
                }
            })
            .finally(() => {
                this.uploadFile = ''
                hideLoader()
            })
        },
        addBuyable: function() {
            showLoader()
            fetch('/api/buyables/add?smiles='+encodeURIComponent(this.addBuyableSmiles)+'&ppg='+this.addBuyablePrice+'&source='+this.addBuyableSource+'&allowOverwrite='+this.allowOverwrite)
            .then(resp => resp.json())
            .then(json => {
                if (json.error || json.success != true) {
                    alert('Error adding buyable compound')
                }
                else {
                    if (json.inserted) {
                        this.buyables.unshift(json.inserted)
                    }
                    if (json.updated) {
                        for (buyable of this.buyables) {
                            if (buyable._id == json.updated._id) {
                                buyable.ppg = json.updated.ppg
                                buyable.source = json.updated.source
                            }
                        }
                    }
                    if (json.duplicate) {
                        alert('Compound already exists in database! Check allow overwrite checkbox to allow overwriting!')
                    }
                }
            })
            .finally(() => hideLoader())
        },
        deleteBuyable: function(_id) {
            if (!window.confirm('Click "OK" to confirm deleting this entry')) {
                return
            }
            showLoader()
            fetch('/api/buyables/delete?_id='+encodeURIComponent(_id))
            .then(resp => resp.json())
            .then(json => {
                if (json.error) {
                    alert(json.error)
                }
                if (json['success']) {
                    for (i in this.buyables) {
                        if (this.buyables[i]['_id']==_id) {
                            this.buyables.splice(i, 1)
                        }
                    }
                }
            })
            .finally(() => hideLoader())
        }
    },
    delimiters: ['%%', '%%'],
});