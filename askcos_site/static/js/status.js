var app = new Vue({
    el: '#app',
    data: {queues: []},
    created: function() {
        this.update();
    },
    methods: {
        update: function() {
            showLoader();
            fetch('/api/v2/celery/')
            .then(resp => resp.json())
            .then(json => {
                console.log(json['queues']);
                this.queues = json['queues'];
            })
            .finally(() => hideLoader())
        }
    },
    delimiters: ['%%', '%%'],
});