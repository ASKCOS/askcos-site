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

new Vue({
    el: "#app",
    data: {
        activeItem: 'chemicals',
        chemicals: [],
        reactions: [],
    },
    created: function() {
        this.loadChemicals();
        this.loadReactions();
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
