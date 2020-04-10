function getCookie(cname) {
    var name = cname + "=";
    var cookie_str = document.cookie;
    if (cookie_str && cookie_str !== '') {
        var cookie_parts = cookie_str.split(';');
        for ( var i = 0; i <cookie_parts.length; i++ ) {
            var c = cookie_parts[i].trim();
            if (c.indexOf(name) === 0) {
                return decodeURIComponent(c.substring(name.length, c.length));
            }
        }
    }
    return undefined;
}

Vue.component('modal', {
    template: '#modal-template'
});

Vue({
    el: "#app",
    data: {
        activeItem: 'chemicals',
        chemicals: [],
        reactions: [],
        showAddModal: false,
        newSmiles: '',
        newDesc: '',
        newActive: true,
        newType: 'chemicals',
        filterActive: 'all',
    },
    created: function() {
        this.loadChemicals();
        this.loadReactions();
    },
    computed: {
        filteredChemicals: function () {
            if (this.filterActive === 'active') {
                return this.chemicals.filter(entry => entry.active)
            } else if (this.filterActive === 'inactive') {
                return this.chemicals.filter(entry => !entry.active)
            } else {
                return this.chemicals
            }
        },
        filteredReactions: function () {
            if (this.filterActive === 'active') {
                return this.reactions.filter(entry => entry.active)
            } else if (this.filterActive === 'inactive') {
                return this.reactions.filter(entry => !entry.active)
            } else {
                return this.reactions
            }
        },
    },
    methods: {
        isActive: function (menuItem) {
            return this.activeItem === menuItem
        },
        setActive: function (menuItem) {
            this.activeItem = menuItem
        },
        loadChemicals: function () {
            this.loadCollection('chemicals')
        },
        loadReactions: function () {
            this.loadCollection('reactions')
        },
        loadCollection: function (category) {
            fetch(`/api/v2/blacklist/${category}/`)
            .then(resp => resp.json())
            .then(json => {
                json.forEach(function (doc) {
                    doc.created = dayjs(doc.created).format('MMMM D, YYYY h:mm A');
                });
                if ( category === 'chemicals') {
                    this.chemicals = json
                } else {
                    this.reactions = json
                }
            })
        },
        addEntry: function () {
            var body = {
                smiles: this.newSmiles,
                description: this.newDesc || 'no description',
                active: this.newActive,
            };
            fetch(`/api/v2/blacklist/${this.newType}/`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': getCookie('csrftoken')
                },
                body: JSON.stringify(body)
            })
                .then(resp => {
                    if (!resp.ok) {
                        throw Error(resp.statusText)
                    }
                    return resp.json();
                })
                .then(json => {
                    json.created = dayjs(json.created).format('MMMM D, YYYY h:mm A');
                    this.showAddModal = false;
                    if ( this.newType === 'chemicals') {
                        this.chemicals.push(json);
                        this.activeItem = 'chemicals'
                    } else {
                        this.reactions.push(json);
                        this.activeItem = 'reactions'
                    }
                })
                .catch(error => console.log(error))
        },
        deleteChemical: function (id) {
            this.deleteEntry(id, 'chemicals')
        },
        deleteReaction: function (id) {
            this.deleteEntry(id, 'reactions')
        },
        deleteEntry: function (id, category) {
            if (window.confirm('Click "OK" to confirm deleting this entry')) {
                fetch(`/api/v2/blacklist/${category}/${encodeURIComponent(id)}`, {
                    method: 'delete',
                    headers: {
                        'X-CSRFToken': getCookie('csrftoken'),
                    },
                })
                    .then(resp => resp.json())
                    .then(json => {
                        if (json['success']) {
                            var collection;
                            if ( category === 'chemicals') {
                                collection = this.chemicals
                            } else {
                                collection = this.reactions
                            }
                            for ( var i = 0; i < collection.length; i++ ) {
                                if (collection[i]['id'] === id) {
                                    collection.splice(i, 1)
                                }
                            }
                        }
                    })
            }
        },
        activateChemical: function (id) {
            this.toggleActivation(id, 'chemicals', 'activate')
        },
        deactivateChemical: function (id) {
            this.toggleActivation(id, 'chemicals', 'deactivate')
        },
        activateReaction: function (id) {
            this.toggleActivation(id, 'reactions', 'activate')
        },
        deactivateReaction: function (id) {
            this.toggleActivation(id, 'reactions', 'deactivate')
        },
        toggleActivation: function (id, category, action) {
            fetch(`/api/v2/blacklist/${category}/${encodeURIComponent(id)}/${action}/`)
            .then(resp => resp.json())
            .then(json => {
                if (json['success']) {
                    var updatedEntry = json['data'];
                    updatedEntry.created = dayjs(updatedEntry.created).format('MMMM D, YYYY h:mm A');
                    var collection;
                    if ( category === 'chemicals') {
                        collection = this.chemicals
                    } else {
                        collection = this.reactions
                    }
                    for ( var i = 0; i < collection.length; i++ ) {
                        if (collection[i]['id'] === id) {
                            collection.splice(i, 1, updatedEntry)
                        }
                    }
                }
            })
        },
    },
    delimiters: ['%%', '%%'],
})
